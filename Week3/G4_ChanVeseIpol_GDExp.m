function [ phi ] = G4_ChanVeseIpol_GDExp( I, phi_0, mu, nu, eta, lambda1, lambda2, tol, epHeaviside, dt, iterMax, reIni )
%Implementation of the Chan-Vese segmentation following the explicit
%gradient descent in the paper of Pascal Getreur "Chan-Vese Segmentation".
%It is the equation 19 from that paper

% I     : Gray color image to segment
% phi_0 : Initial phi
% mu    : mu lenght parameter (regularizer term)
% nu    : nu area parameter (regularizer term)
% eta   : epsilon for the total variation regularization
% lambda1, lambda2: data fidelity parameters
% tol   : tolerance for the sopping criterium
% epHeaviside: epsilon for the regularized heaviside.
% dt     : time step
% iterMax : Maximum number of iterations
% reIni   : Iterations for reinitialization. 0 means no reinitialization

[ni,nj] = size(I);
% Definition of index vectors for the inner points
i = 2:ni+1;
j = 2:nj+1;
% step size for difference functions
hi=1;
hj=1;

% initializations for the loop. Space for the ghost boundaries is
% allocated.
phi(ni+2,nj+2) = 0;
phi(i,j) = phi_0;
dif = inf;
nIter = 0;
% The operation is terminated if the difference between the consecutive
% level set functions are below a certain value or a number of iterations
% are completed
while dif>tol && nIter<iterMax
  % The loop variables are updated
  phi_old = phi;
  nIter = nIter+1;
  phi_in = phi(2:end-1,2:end-1);
  %Fixed phi, Minimization w.r.t c1 and c2 (constant estimation)
  % The values of the regions inside and outside of the curve are updated
  % using the average intensities of the image. As the regularization is
  % not necessary for this calculation, the discretisized versions
  % continuous of the expressions are not preferred.
  c1 = sum(sum(I(phi_in>=0))) / sum(sum((phi_in>=0)));
  c2 = sum(sum(I(phi_in<0))) / sum(sum((phi_in<0)));

  %Boundary conditions (ghost boundaries are added for difference
  % operations)
  phi(1,:)   = phi_old(2,:);
  phi(end,:) = phi_old(end-1,:);

  phi(:,1)   = phi_old(:,2);
  phi(:,end) = phi_old(:,end-1);

  %Regularized Dirac's Delta computation
  delta_phi = G4_diracReg(phi, epHeaviside);   %notice delta_phi=H'(phi)

  %derivatives estimation
  %i direction, forward finite differences
  phi_iFwd  = DiFwd(phi_old, hi);
  phi_iBwd  = DiBwd(phi_old, hi);

  %j direction, forward finite differences
  phi_jFwd  = DjFwd(phi_old, hj);
  phi_jBwd  = DjBwd(phi_old, hj);

  %centered finite differences
  phi_icent = (phi_iFwd+phi_iBwd)/2;
  phi_jcent = (phi_jFwd+phi_jBwd)/2;

  %A and B estimation (from the Pascal Getreuer's IPOL paper "Chan Vese
  % Segmentation"
  A = mu./sqrt(eta^2 + phi_iFwd.^2 + phi_jcent.^2);
  B = mu./sqrt(eta^2 + phi_icent.^2 + phi_jFwd.^2);

  %Equation 22. The image indices are decreased as adding ghost boundaries
  % to the image is not necessary.
  phi(i, j) = ...
    (phi_old(i, j) + dt .* delta_phi(i, j) ...
    .*(A(i, j) .* phi_old(i+1, j)  + A(i-1, j) .* phi(i-1, j) ...
    +  B(i, j) .* phi_old(i, j+1)  + B(i, j-1) .* phi(i, j-1) ...
    - nu - lambda1.*((I(i-1, j-1)-c1).^2) +...
           lambda2.*((I(i-1, j-1)-c2).^2))) ...
    ./(1 + dt.*delta_phi(i, j).*...
    (A(i, j) + A(i-1, j) + B(i, j) + B(i, j-1)));

  %Reinitialization of phi to ensure stable evolution of the curve(s)
  if reIni>0 && mod(nIter, reIni)==0
    indGT = phi >= 0;
    indLT = phi <  0;

    phi=double(bwdist(indLT) - bwdist(indGT));

    %Normalization [-1 1]
    nor = min(abs(min(phi(:))), max(phi(:)));
    phi=phi/nor;
  end

  %Difference. This stopping criterium has the problem that phi can
  %change but not the zero level set, that it really is what we are
  %looking for.
  dif = mean(sum((phi(:) - phi_old(:)).^2));

  %Plot the level sets surface
  subplot(1,2,1)
  %The level set function
  surfc(phi)
  hold on
  %The zero level set over the surface
  contour(phi>0);
  hold off
  title('Phi Function');

  %Plot the curve evolution over the image
  subplot(1,2,2)
  imagesc(I);
  colormap gray;
  hold on;
  contour(phi>0,'r')
  title('Image and zero level set of Phi')

  axis off;
  hold off
  drawnow;
  pause(.0001);
end
% To delete the ghost boundaries
phi=phi(i,j);