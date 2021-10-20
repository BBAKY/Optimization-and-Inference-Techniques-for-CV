function [u] = G4_Poisson_Equation_Axb(f, dom2Inp, param)
%this code is not intended to be efficient. 
%f: original image
%dom2Inp: mask specifying the pixels to be inpainted
[ni, nj]=size(f);

%We add the ghost boundaries (for the boundary conditions)
f_ext = zeros(ni+2, nj+2);
f_ext(2:end-1, 2:end-1) = f;
dom2Inp_ext =zeros(ni+2, nj+2);
dom2Inp_ext (2:end-1, 2:end-1) = dom2Inp;

%Store memory for the A matrix and the b vector    
nPixels =(ni+2)*(nj+2); %Number of pixels together with ghost boundaries

%We will create A sparse, this is the number of nonzero positions

%idx_Ai: Vector for the nonZero i index of matrix A, corresponds to the
%equation number
%idx_Aj: Vector for the nonZero j index of matrix A, corresponds to pixels
%of the new image
%a_ij: Vector for the value at position ij of matrix A, corresponds to the
%coefficients of the pixel values of the new images

%Preallocation of the RHS vector which will hold the original pixel values
b = zeros(nPixels,1);

%Vector counter for sparse coefficient matrix
%It is increased for each coefficient definiton
idx=1;

%Ghost boundaries are equated to the boundaries of the original image as a
%part of the solution in case pixels on the boundary are to be inpainted
%Corner pixel equations, which are not useful, for meeting boundary sides 
%combine to make an equation for each corner. Therefore, the
%coefficient matrix size is not affected. 

%North side boundary conditions
i=1;
for j=1:nj+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate,
    %increases column by column
    p = (j-1)*(ni+2)+i;

    idx_Ai(idx) = p;
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx=idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p+1;% i increased for the lower neighbour
    a_ij(idx) = -1;   
    idx=idx+1;
    %There is nothing independent of the pixel values in the equation, 
    %added for completeness as they are already zero due to preallocation        
    b(p) = 0;
end

%South side boundary conditions
i=ni+2;
for j=1:nj+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;
    
    idx_Ai(idx)= p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx=idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p-1;% i decreased for the upper neighbour
    a_ij(idx) = -1;   
    idx=idx+1;
            
    b(p) = 0;
    
end

%West side boundary conditions
j=1;
for i=1:ni+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;

    idx_Ai(idx)=p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx=idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p+(ni+2);%j increased for the right neighbour
    a_ij(idx) = -1;   
    idx=idx+1;
            
    b(p) = 0;
    
    
end

%East side boundary conditions
j=nj+2;
for i=1:ni+2
    %from image matrix (i,j) coordinates to vectorial (p) coordinate
    p = (j-1)*(ni+2)+i;
    
    idx_Ai(idx)=p; 
    idx_Aj(idx) = p; 
    a_ij(idx) = 1;
    idx=idx+1;
    
    idx_Ai(idx) = p;
    idx_Aj(idx) = p-(ni+2);% j decreased for the left neighbour
    a_ij(idx) = -1;   
    idx=idx+1;
            
    b(p) = 0;
    
end

%Inner points are equated to the original or computed from their neighbours
%As the solution is implicit, pixels are not necessarily isolated
for j=2:nj+1
    for i=2:ni+1
     
        %from image matrix (i,j) coordinates to vectorial (p) coordinate
        p = (j-1)*(ni+2)+i;
                                            
        if (dom2Inp_ext(i,j)==1) %If we have to inpaint this pixel
            
            idx_Ai(idx)=p; 
            idx_Aj(idx) = p-(ni+2); 
            a_ij(idx) = -1;
            idx=idx+1;
            
            idx_Ai(idx) = p;
            idx_Aj(idx) = p-1;
            a_ij(idx) = -1;   
            idx=idx+1;
            % for the center pixel 
            idx_Ai(idx)=p; 
            idx_Aj(idx) = p; 
            a_ij(idx) = 4;
            idx=idx+1;
            
            idx_Ai(idx) = p;
            idx_Aj(idx) = p+1;
            a_ij(idx) = -1;   
            idx=idx+1;
            
            idx_Ai(idx) = p;
            idx_Aj(idx) = p+(ni+2);
            a_ij(idx) = -1;   
            idx=idx+1;
                    
            if (isfield(param, 'driving'))
                b(p) = param.driving(i-1,j-1);
            else
                b(p) = 0;
            end
    
        else %we do not have to inpaint this pixel 
            
            idx_Ai(idx) = p;
            idx_Aj(idx) = p;
            a_ij(idx) = 1;   
            idx=idx+1;
            %The values  of the pixels not specified by the mask (to copy 
            %from the original image)
            b(p) = f_ext(i,j);
             
        end       
    end
end
    %A is a sparse matrix, so for memory requirements we create a sparse
    %matrix

    A = sparse(idx_Ai, idx_Aj, a_ij, nPixels, nPixels);
    
    %Solve the system of equations
    x = mldivide(A,b);
    
    %From vector to matrix (of the image)
    u_ext = reshape(x, ni+2, nj+2);
    
    %Eliminate the ghost boundaries as they are not the part of the image
    u = full(u_ext(2:end-1, 2:end-1));