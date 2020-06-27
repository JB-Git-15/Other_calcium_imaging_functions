//// Compute dissimilarity matrix 
////  1 - before using do (only the first time):  mex  Dissimilarity_mex.c
////  2 - then call it from matlab:  Dissimilarity = Dissimilarity_mex(Matrix)
////      Matrix is a N cell by N stimuli by N repetitions
////      JB 26/06/2020


#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        mxArray *Trace_in, *Dissimilarity_out, *Rho_temp, *Ind_stims_tempU, *Ind_stims_tempD, *Ind_stims_Diag_temp, *PairsOfStims_temp;
        const mwSize *dims;
        double *Trace, *Dissimilarity, *Rho, *Ind_stimsU, *Ind_stimsD,*Ind_stims_Diag, *PairsOfStims; 
        int N_cells, N_stims, N_rep;
        double m1,m2, std1, std2, cov, c ; // means and standard deviations
        int n_stim1, n_stim2, n_stim3;
                
        // Input array:
        Trace_in = mxDuplicateArray(prhs[0]);

        // Get dimensions of input array:
        dims    = mxGetDimensions(prhs[0]);
        N_cells = (int) dims[0];
        N_stims = (int) dims[1];
        N_rep   = (int) dims[2];
        
        // Output array:
        Dissimilarity_out   = plhs[0] = mxCreateDoubleMatrix(N_stims, N_stims, mxREAL);
        Rho_temp            = mxCreateDoubleMatrix(N_stims, N_stims, mxREAL);
        Ind_stims_tempU     = mxCreateDoubleMatrix(N_stims*(N_stims-1)/2, 1, mxREAL);
        Ind_stims_tempD     = mxCreateDoubleMatrix(N_stims*(N_stims-1)/2, 1, mxREAL);
        Ind_stims_Diag_temp = mxCreateDoubleMatrix(N_stims,2,mxREAL);
        PairsOfStims_temp   = mxCreateDoubleMatrix(N_stims,2,mxREAL);
        
        // Access the contents of the input and output arrays:
        Trace               = mxGetPr(Trace_in);
        Rho                 = mxGetPr(Rho_temp);
        Dissimilarity       = mxGetPr(Dissimilarity_out);
        Ind_stimsU          = mxGetPr(Ind_stims_tempU);
        Ind_stimsD          = mxGetPr(Ind_stims_tempD);
        Ind_stims_Diag      = mxGetPr(Ind_stims_Diag_temp);
        PairsOfStims        = mxGetPr(PairsOfStims_temp);
        
        // list the indices of the triangular superior elements of the N_stims by N_stims matrix
        // also get the indexes of the diagonal elements of the N_stims by N_stims matrix
        int count1 = 0; 
        int count2 = 0;
        for (int col = 0; col < N_stims; col++) {
            for (int lign = 0; lign < N_stims; lign++) {
                if (col >lign){
                      Ind_stimsU[count1]             = (double)(col*N_stims  + lign);
                      Ind_stimsD[count1]             = (double)(lign*N_stims + col);
                      PairsOfStims[count1]           = (double) lign;
                      PairsOfStims[count1 + N_stims] = (double) col;
                      count1 = count1 + 1;
                }
                else if (col == lign) {
                      Ind_stims_Diag[count2]         = (double)(lign*N_stims + col);
                      count2                         =  count2 +1;
                }                
            }
        }       
        
        for (int n_pair_stims = 0; n_pair_stims < N_stims*(N_stims-1)/2; n_pair_stims++) {
             n_stim1 = (int)PairsOfStims[n_pair_stims];
             n_stim2 = (int)PairsOfStims[n_pair_stims + N_stims];
             c       = 0;
                  for (int n_rep1 = 0; n_rep1 < N_rep; n_rep1++) {
                     for (int n_rep2 = 0; n_rep2 < N_rep; n_rep2++){
                         // compute correlation  :
                              m1 = m2 = std1 = std2 = cov = 0;
                              for (int n_cell = 0; n_cell < N_cells; n_cell++) {
                                 m1 = m1 + Trace[(n_cell + n_stim1*N_cells) + n_rep1*N_stims*N_cells];
                                 m2 = m2 + Trace[(n_cell + n_stim2*N_cells) + n_rep2*N_stims*N_cells];
                               }
                                 m1   = m1/((double)N_cells);
                                 m2   = m2/((double)N_cells);
                              for (int n_cell2 = 0; n_cell2 < N_cells; n_cell2++) {
                                 // printf("Indices a consider %d \n",(int)((n_cell2 + n_stim1*N_cells)+ n_rep1*N_stims*N_cells)); 
                                  std1 = std1 + pow((Trace[(int)((n_cell2 + n_stim1*N_cells)+ n_rep1*N_stims*N_cells)] - m1),2.0);
                                  std2 = std2 + pow((Trace[(int)((n_cell2 + n_stim2*N_cells)+ n_rep2*N_stims*N_cells)] - m2),2.0);
                                  cov  = cov  + (Trace[(int)((n_cell2 + n_stim1*N_cells)+ n_rep1*N_stims*N_cells)] - m1)*
                                                (Trace[(int)((n_cell2 + n_stim2*N_cells)+ n_rep2*N_stims*N_cells)] - m2);
                              } 
                                 std1 = sqrt(std1/((double)N_cells));
                                 std2 = sqrt(std2/((double)N_cells));
                                 cov  = cov/((double)N_cells);  
                         c = c +  cov/(std1*std2);  
                     }
                  }   
           c       =  c/(N_rep*N_rep);
           Rho[(int)Ind_stimsU[n_pair_stims]]  = c;
           Rho[(int)Ind_stimsD[n_pair_stims]]  = c;  // we directly transpose 
        }       
                
        for (int n_pair_stims1 = 0; n_pair_stims1 < N_stims; n_pair_stims1++) {
           c      = 0;
             for (int n_rep3 = 0; n_rep3 < N_rep; n_rep3++) {
                for (int n_rep4 = 0; n_rep4 < N_rep; n_rep4++){
                    if (n_rep3 != n_rep4){
                       m1 = m2 = std1 = std2 = cov = 0;
                     // compute correlation
                       for (int n_cell3 = 0; n_cell3 < N_cells; n_cell3++) {
                                 m1 = m1 + Trace[(n_cell3 + n_pair_stims1*N_cells) + n_rep3*N_stims*N_cells];
                                 m2 = m2 + Trace[(n_cell3 + n_pair_stims1*N_cells) + n_rep4*N_stims*N_cells];
                               }
                                 m1   = m1/((double)N_cells);
                                 m2   = m2/((double)N_cells);
                       for (int n_cell4 = 0; n_cell4 < N_cells; n_cell4++) {
                                 // printf("Indices a consider %d \n",(int)((n_cell2 + n_stim1*N_cells)+ n_rep1*N_stims*N_cells)); 
                                  std1 = std1 + pow((Trace[(int)((n_cell4 + n_pair_stims1*N_cells)+ n_rep3*N_stims*N_cells)] - m1),2.0);
                                  std2 = std2 + pow((Trace[(int)((n_cell4 + n_pair_stims1*N_cells)+ n_rep4*N_stims*N_cells)] - m2),2.0);
                                  cov  = cov  + (Trace[(int)((n_cell4 + n_pair_stims1*N_cells)+ n_rep3*N_stims*N_cells)] - m1)*
                                                (Trace[(int)((n_cell4 + n_pair_stims1*N_cells)+ n_rep4*N_stims*N_cells)] - m2);
                              } 
                                 std1 = sqrt(std1/((double)N_cells));
                                 std2 = sqrt(std2/((double)N_cells));
                                 cov  = cov/((double)N_cells);                  
                                 c    = c + cov/(std1*std2); 
                    }    
                }
             } 
            Rho[(int)Ind_stims_Diag[n_pair_stims1]]  = c/(N_rep*(N_rep-1));  // diagonal term  
        }
        
         /*   for (int ligne_1 = 0; ligne_1 < N_stims; ligne_1++) {
                 for (int col_1 = 0; col_1 < N_stims; col_1++) {  
                      if (ligne_1 != col_1) {
                          Dissimilarity[(int)(ligne_1 + col_1*N_stims)]  
          *                                 =  .5*(1 - Rho[(int)(ligne_1 + col_1*N_stims)]/sqrt(abs(Rho[(int)(ligne_1 + ligne_1*N_stims)]
                                                                                   *Rho[(int)(col_1 + col_1*N_stims)]))); 
                      } 
                      else{
                          Dissimilarity[(int)(ligne_1 + ligne_1*N_stims)] = 0;
                      }
                 }
               }
          */
        
        for (int ligne_x = 0; ligne_x < N_stims*N_stims; ligne_x++) {
           Dissimilarity[ligne_x] = Rho[ligne_x];
        
        }
        
  }
        
        
        
        
        
        
                
                 
                
          
                                         
                         

              
                       

           
       
        
        
        
        
        
/* In matlab  : input : Trace 
       Index_d              = find(diag(ones(Nstims,1)));
       for n_stims = 1:Nstims
               c = 0;    
           for n_rep_1 = 1:Reps
               Int          = 1:Reps;
               Int(n_rep_1) = [];
               for n_rep_2  = 1:length(Int)
                   n_rep_3  = Int(n_rep_2);                  
                   c        = c + corr(Trace(:,n_stims,n_rep_1),Trace(:,n_stims,n_rep_3));
               end    
           end
                   c        = c/(Reps*(Reps-1)); 
                   Rho(Index_d(n_stims)) = c;            
       end
       
        
       
           for n_1 = 1:Nstims 
                 Int = 1:Nstims;
                 Int(n_1) = [];
                for n_2 = 1:length(Int)
                     n_3 = Int(n_2);
                     Dissimilarity(n_1,n_3) = .5*(1 - Rho(n_1,n_3)/sqrt(abs(Rho(n_1,n_1)*Rho(n_3,n_3))));
                end    
           end

*/






