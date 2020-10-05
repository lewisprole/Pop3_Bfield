int merged[NSinksAllTasks] = { 0 };  /* keep track of deleted sinks */
struct merged_SinkP;   /* will overwrite SinkP struct after mergers */
  {
    double Mass;
    double Posx;
    double Posy;
    double Posz;
    double Velx;
    double Vely;
    double Velz;
  }
int N_survived = 0 ;  /*counter */

for(j = 0; j < NSinksAllTasks; j++)
  {
    if (merged[j]==1)  /* if it no longer exists, do nothing */
      break;
    for(i = 0; i < NSinksAllTasks; i++)
      {
        int to_be_merged[NSinksAllTasks] = { 0 };  /* tracks which sinks are to be merged with sink [j]*/
        if (merged[i]==1)
          break;  /* can’t merge with a sink that doesn’t exist anymore */

        if (SinkP[i].ID ==  SinkP[j].ID)
          break;  /* can’t merge with itself */
       
        /* Check separation*/
        dx   = GRAVITY_NEAREST_X(SinkP[i].Pos[0] - SinkP [j].Pos[0]);
        dy   = GRAVITY_NEAREST_Y(SinkP[i].Pos[1] - SinkP [j].Pos[1]);
        dz   = GRAVITY_NEAREST_Z(SinkP[i].Pos[2] - SinkP [j].Pos[2]);
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

         ekin = 0.5 * P[icell].Mass * dv; /*energies*/
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
         else 
           to_be_merged [i]=1; /*passed tests - good to merge*/

       } /*end of [i] loop*/

     /*Perform merge - properties of sink [i] transfered to sink [j]*/
     for (int i = 0; i < NSinksAllTasks; i++)                
       {
         if (to_be_merged[i]==1)
           {
             /*move to centre of mass */
             SinkP[j].Pos[0]=( SinkP[j].Mass* SinkP[j].Pos[0] + SinkP[i].Mass* SinkP[i].Pos[0] ) / (SinkP[j].Mass+ SinkP[i].Mass);
             SinkP[j].Pos[1]=( SinkP[j].Mass* SinkP[j].Pos[1] + SinkP[i].Mass* SinkP[i].Pos[1] ) / (SinkP[j].Mass+ SinkP[i].Mass);
             SinkP[j].Pos[2]=( SinkP[j].Mass* SinkP[j].Pos[2] + SinkP[i].Mass* SinkP[i].Pos[2] ) / (SinkP[j].Mass+ SinkP[i].Mass);

             /*transfer linear momentum*/
             SinkP[j]. Vel[0] = (SinkP[j].Mass * SinkP[j]. Vel[0] + SinkP[i].Mass * SinkP[i]. Vel[0]) / (SinkP[j].Mass+ SinkP[i].Mass);
             SinkP[j]. Vel[1] = (SinkP[j].Mass * SinkP[j]. Vel[1] + SinkP[i].Mass * SinkP[i]. Vel[1]) / (SinkP[j].Mass+ SinkP[i].Mass);
             SinkP[j]. Vel[0] = (SinkP[j].Mass * SinkP[j]. Vel[2] + SinkP[i].Mass * SinkP[i]. Vel[2]) / (SinkP[j].Mass+ SinkP[i].Mass);

             /*transfer mass*/
             SinkP[j].Mass += SinkP[i].Mass;
             merged[i]=1; /*let the [j] loop know which sinks no longer exist*/
           {
       {
     }     /*end of the [j] loop */     

  /*store surviving sink data, get rid of merged sinks by storing data in a different struct containing only surviving sinks */     
  for (int i = 0; i < NSinksAllTasks; i++)
    {
       if merged[i]=0; /*survived*/
       merged_SinkP[N_survived].Mass = SinkP[i].Mass;
       merged_SinkP[N_survived].Posx = SinkP[i].Pos[0];
       merged_SinkP[N_survived].Posy = SinkP[i].Pos[1];
       merged_SinkP[N_survived].Posz = SinkP[i].Pos[2];
       merged_SinkP[N_survived].Velx = SinkP[i].Vel[0];
       merged_SinkP[N_survived].Vely = SinkP[i].Vel[1];
       merged_SinkP[N_survived].Velz = SinkP[i].Vel[2];
       N_survived += 1;
     }  


