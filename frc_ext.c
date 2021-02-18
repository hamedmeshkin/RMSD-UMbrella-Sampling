/* ##   The CA of the resids that participate in the bias RMSD potential are:
##### 360 to 370. This script first align the namd simulated coordiante of CAs of resid 353 354 360 to 370 with 
##### initial coordinate of these CAs at 2JLN and 4D1B sructure, then RMSD of CAs of resid 360 to 370 measue by:
#####
##### Reaction Coordinate = (RMSD_2JLN / RMS + RMSD_4D1B / RMS) + 0.5  
##### The Reaction Coordinate is a vaue from 0.0 to 1.0 which 0.0 is 2JLN and 1.0 is 4D1B state
##### RMSD_2JLN is alignment of the simulated coordinate with 2JLN structure
##### RMSD_4D1B is alignment of the simulated coordinate with 4D1B structure
##### RMS is the RMSD of 2JLN and 4D1B
## */
 #include <math.h>
 #include <iostream>
 #include <tcl.h>
 #include <stdlib.h>
 #include <stdio.h>
 #include <complex.h>
 #include "Eigen/Eigen"
 
 
 static double  Const, RmRef, RMSD4, RMSD2, RMSD, KwRM, RMS, Coeff, DeltaRmsd; 
 static double  coor4[13][3], coor2[13][3], MR2[3][13] = {0} , MR4[3][13] = {0};   
 static float   stop = 1;  
 static char    buffer1 [20], buffer3 [100], buffer5 [10];
 static int     sysiz, Rindex[13], replica_id, i_job, ts, oldTS = -1;
 int Tables();
 static FILE    *gFile; 
  
 using namespace Eigen;
 using namespace std;
 
// ##############################################################################################################
// ############################### Finding Rotation and Translation Matrix by SVD################################
// ##############################################################################################################
 Affine3d Find3DAffineTransform(Matrix3Xd P, Matrix3Xd Q) {

    // Default output
    Affine3d A;
    A.linear() = Matrix3d::Identity(3, 3);
    A.translation() = Vector3d::Zero();

    if (P.cols() != Q.cols())
        throw "Find3DAffineTransform(): input data mis-match";

    // Center the data
    Vector3d p = P.rowwise().mean();
    Vector3d q = Q.rowwise().mean();

    Matrix3Xd X = P.colwise() - p;
    Matrix3Xd Y = Q.colwise() - q;

    // SVD
    MatrixXd Cov = X*Y.transpose();
    JacobiSVD<MatrixXd> svd(Cov, ComputeThinU | ComputeThinV);

    // Find the rotation, and prevent reflections
    Matrix3d IiI = Matrix3d::Identity(3, 3);
    double d = (svd.matrixV()*svd.matrixU().transpose()).determinant();
    (d > 0.0) ? d = 1.0 : d = -1.0;
    IiI(2, 2) = d;
  
    Matrix3d R = svd.matrixV()*IiI*svd.matrixU().transpose();

    // The final transform
    A.linear() = R;
    A.translation() = q - R*p;

    return A;
}

// ##############################################################################################################
// ############################### Main Script for Hormonic Potential on RMSD ###################################
// ##############################################################################################################
 extern "C++" int rms_CNT(ClientData clientData, Tcl_Interp *interp, 
        int objc, Tcl_Obj *CONST objv[]) {
     
    static Matrix <double,3,13> oldpoints4,oldpoints2,newpoints,remappedpoints4,remappedpoints2,r4,r2;   
    double  X,Y,Z ;
    int     num=0, i;
    Tcl_Obj   **list;
    Tcl_Obj   *Coord1;
        
    Tcl_Obj *value1 = Tcl_GetVar2Ex(interp, "RmRef", NULL, TCL_GLOBAL_ONLY);
    Tcl_GetDoubleFromObj(interp,value1,&RmRef); // RMSD Refrence Point     
    Tcl_Obj *value3 = Tcl_GetVar2Ex(interp, "ts", NULL, TCL_GLOBAL_ONLY);
    Tcl_GetIntFromObj(interp,value3,&ts);
    
    // Getting the initial value of the script which are read only once    
    if (stop == 1) {

        stop = 0;
        
        Tcl_Obj *value2 = Tcl_GetVar2Ex(interp, "replica_id", NULL, TCL_GLOBAL_ONLY);
        Tcl_GetIntFromObj(interp,value2,&replica_id);
        Tcl_Obj *value4 = Tcl_GetVar2Ex(interp, "i_job", NULL, TCL_GLOBAL_ONLY);
        Tcl_GetIntFromObj(interp,value4,&i_job);
        Tcl_Obj *value6 = Tcl_GetVar2Ex(interp, "RMS", NULL, TCL_GLOBAL_ONLY);
        Tcl_GetDoubleFromObj(interp,value6,&RMS);  // Rmsd of 2JLN and 4D1B
        Tcl_Obj *value5 = Tcl_GetVar2Ex(interp, "KwRM", NULL, TCL_GLOBAL_ONLY);
        Tcl_GetDoubleFromObj(interp,value5,&KwRM); // Spring Constant
        Tcl_Obj *value3 = Tcl_GetVar2Ex(interp, "ts", NULL, TCL_GLOBAL_ONLY);
        Tcl_GetIntFromObj(interp,value3,&ts);
        
        Tables();
        
        for(unsigned int i = 0; i < sysiz; ++i) {
            Vector3d Points_One( coor4[i][0], coor4[i][1], coor4[i][2]);
            Vector3d Points_Two( coor2[i][0], coor2[i][1], coor2[i][2]);
            oldpoints4.col(i)=Points_One;
            oldpoints2.col(i)=Points_Two;
        }
        cout << "oldpoints4 =" << endl << oldpoints4 << endl;
        cout << "oldpoints2 =" << endl << oldpoints2 << endl;
        RMS = RMS * sqrt(sysiz-2.0) * 2.0;
        Coeff = KwRM / RMS; 
    }   

    // Getting the coordinate of the atoms at each dt of simulation
    for (i=0 ; i < sysiz ; i++) { 
        snprintf ( buffer1, 20, "atmcrd(%d)" , Rindex[i] );        
        Coord1 = Tcl_GetVar2Ex(interp, buffer1 , NULL, TCL_GLOBAL_ONLY);
        if (Tcl_ListObjGetElements(interp, Coord1, &num, &list) != TCL_OK) return TCL_ERROR;
        
        if (Tcl_GetDoubleFromObj(interp, list[0], &X) != TCL_OK) return TCL_ERROR;
        if (Tcl_GetDoubleFromObj(interp, list[1], &Y) != TCL_OK) return TCL_ERROR;
        if (Tcl_GetDoubleFromObj(interp, list[2], &Z) != TCL_OK) return TCL_ERROR;
        Vector3d Points_Sim( X, Y, Z);
        newpoints.col(i)=Points_Sim;
    }

    // Use SVD alignment algorithm
    Affine3d A4 = Find3DAffineTransform(oldpoints4,newpoints);
    Affine3d A2 = Find3DAffineTransform(oldpoints2,newpoints);
    
    // apply transformation to restore the points
    remappedpoints4 = (A4.linear()*oldpoints4).colwise() + A4.translation();
    remappedpoints2 = (A2.linear()*oldpoints2).colwise() + A2.translation();
    
    // Measuring the RMSD
    r4 = newpoints - remappedpoints4;
    r2 = newpoints - remappedpoints2;
    
    double rmsd4 = 0, rmsd2 = 0;
    for(unsigned int i = 0; i < 3; ++i) { // i has 1,2 and 3 value which are representing as x y z  
        for(unsigned int j = 2; j < sysiz; ++j) { // j starts from 2 to ignore CA of 353 and 354
            MR4[i][j] = r4(i,j);
            MR2[i][j] = r2(i,j);
            rmsd4 += (MR4[i][j] * MR4[i][j]);
            rmsd2 += (MR2[i][j] * MR2[i][j]);
        }
    }
    
    RMSD4 = sqrt(rmsd4);
    RMSD2 = sqrt(rmsd2);
    RMSD = 0.5 + (RMSD2-RMSD4) / RMS;
    
    // Writing to a File and pass the variables to Tcl
    if (oldTS < ts){
        fprintf (gFile,"%0.5f\n",RMSD);
        oldTS = ts;
    }
    
    if (ts % 50000 == 0) {
        fflush (gFile);
        if (fflush(gFile)) { 
          perror ("Flush Error");
          exit (EXIT_FAILURE);
        }
    }
    
    ts++;
    snprintf( buffer5, 10, "%f" ,RMSD);
    Tcl_SetVar(interp, "RMSD", buffer5, TCL_GLOBAL_ONLY);
    
    return TCL_OK;    
}

// #########################################################################################################
// ############################### Applying RMSD harmonic Force on Atoms ###################################
// #########################################################################################################
extern "C++" int rmsF_CNT(ClientData clientData, Tcl_Interp *interp, 
        int objc, Tcl_Obj *CONST objv[]) {
    
    Tcl_Obj *value7 = Tcl_GetVar2Ex(interp, "DeltaRmsd", NULL, TCL_GLOBAL_ONLY);
    Tcl_GetDoubleFromObj(interp,value7,&DeltaRmsd);
 
    Const = Coeff * DeltaRmsd; 
    for (int i=2 ; i < sysiz ; i++) {  // j starts from 2 to ignore CA of 353 and 354
        snprintf ( buffer3, 100, "addforce %d {%.15f %.15f %.15f} " ,
            Rindex[i], Const * (MR4[0][i]/RMSD4 - MR2[0][i]/RMSD2), Const * (MR4[1][i]/RMSD4 - MR2[1][i]/RMSD2), Const * (MR4[2][i]/RMSD4 - MR2[2][i]/RMSD2));
        Tcl_GlobalEval(interp, buffer3);
    }
  
    return TCL_OK;    
}

// #########################################################################################################
// ############################### Assigning a Tcl Name and File For the Script ############################
// #########################################################################################################
extern "C" int Rms_ext_Init(Tcl_Interp *interp) {
    Tcl_CreateObjCommand(interp, "rmsD", rms_CNT,
             (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
    return TCL_OK;
}
extern "C" int Rmsf_ext_Init(Tcl_Interp *interp) {
    Tcl_CreateObjCommand(interp, "rmsDF", rmsF_CNT,
             (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
    return TCL_OK;
}

                     // \\
                    // * \\ 
                   // *** \\ 
                  // ***** \\ 
                 // ******* \\ 
                // ********* \\ 
               // *********** \\ 
              // ************* \\ 
             // *** Tables **** \\ 
            // ***************** \\ 
           // ******************* \\ 
          // ********************* \\ 
         // *********************** \\   
  // *************************************** \\     
 // ****   Index of the Selected Atoms   **** \\   
// ******************************************* \\ 

 static char    buffer4 [35];    
 int Tables() {
    
    FILE *pFile;
    pFile = fopen("input/Rindex.dat","r");
    int i = 0; int num9;
    while (!feof(pFile)) {
        fscanf(pFile,"%d",&num9);
        Rindex[i] = num9;
        i++;
    }
    sysiz = i;
    fclose(pFile);
 
// ******************************************* \\
// ** Coordinate of the selected 2JLN Atoms ** \\
// ******************************************* \\
    
    pFile = fopen("input/coor_2JLN.dat","r");
    i = 0; float num3; int k = 0;
    while (!feof(pFile)) {
        fscanf(pFile,"%fl",&num3);
        if (k == 0) {
            k = 1;
            coor2[i][0] = num3;
        } else if (k == 1) {
            coor2[i][1] = num3;
            k = 2;
        } else  {
            coor2[i][2] = num3;
            k = 0;
            i++;
        }
    }
    fclose(pFile);
 
// ******************************************* \\
// ** Coordinate of the selected 4D1B Atoms ** \\
// ******************************************* \\
    
    pFile = fopen("input/coor_4D1B.dat","r");
    i = 0; float num4;  k = 0;
    while (!feof(pFile)) {
        fscanf(pFile,"%fl",&num4);
        if (k == 0) {
            k = 1;
            coor4[i][0] = num4;
        } else if (k == 1) {
            coor4[i][1] = num4;
            k = 2;
        } else  {
            coor4[i][2] = num4;
            k = 0;
            i++;
        }
    }
    fclose(pFile);   
// ******************************************* \\
// ************** RMSD.dat File ************** \\
// ******************************************* \\
    
    
    snprintf ( buffer4, 35, "output_%03d/dat/RMSD.%03d.dat" ,i_job , replica_id );
    gFile = fopen (buffer4,"w");
 }
                        /*
  // *************************************** \\     
 // ****************** End ****************** \\   
// ******************************************* \\
          \\ *********************** //
           \\ ********************* // 
            \\ ******************* // 
             \\ ***************** // 
              \\ *************** //
               \\ ************* //
                \\ *********** //
                 \\ ********* //
                  \\ ******* // 
                   \\ ***** //
                    \\ *** //
                        */