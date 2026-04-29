function DD=penalty2(M,N,ds) %will work for M not the same as N

%MxN size of image; example: MxN = 64x64, ds=4
%M/ds x N/ds second level = 16x16 (256)
%n0=width of image over which beta parameters are defined (columns)
n0 = N / ds; %16
n0h = M / ds;
n = M*N/ds/ds; %n0h*n0;

dim11 = n0h*(2*n0-1); % no of rows in first matrix (496)
dim12 = 2*n - 2* n0; % no of rows in second matrix (480)

dim1 = dim11+dim12;% no of rows  976 (496 + 480)
dim2 = 2*n; %no of columns (512 for 64x64 images)

DD=sparse(dim1,dim2); %zeros(dim1,dim2);

for i=1:n0h, %build matrix D2
    
    maxim = i*(2*n0-1);
    minim = (i-1)*(2*n0-1) + 1; %maxim - (2*n0-2);
     
    for j=minim:maxim,
        DD(j,j+i-1) = -1;  
        if j+i+1 <= dim2,
            DD(j, j + i+1) = 1;
        end
    end
end


for i=1:dim12, %build matrix D2n0
    DD(i+dim11,i) = -1;
    if i + 2*n0 <=dim2,
        DD(i+dim11, i + 2*n0) = 1;
    end
end
 
