classdef PoissSolv2D < handle
     properties
         MlBas = []; 
         f = [];
         iLvl = [];
         iBasisFctInd = [];
     end
     
      methods (Access = public)
         function obj = PoissSolv2D(MlBas,f)
             if nargin >0
                obj.MlBas = MlBas;
                obj.f = f;
            end
         end
         
         function [Stiffn, rhs, iLvl,iBasisFctIndU,iBasisFctIndV] = assembleMl(objU,objV)
             allKnotsU = objU.MlBas.getAllKnots; % Get knotvector in U direction.
             allKnotsV = objV.MlBas.getAllKnots; % Get knotvector in V direction.
             nOFU = objU.MlBas.nOF; % Number of basis functions in U direction.
             nOFV = objV.MlBas.nOF; % Number of basis functions in V direction.
             nOF = nOFU*nOFV; % Number of basis function in the tensor product U*V.
             nOEU =length(allKnotsU) -1; % number of elements
             nOEV =length(allKnotsV) -1; % number of elements
 
             % Find basis functions.
%              Nx = objU.MlBas.levelBas{1}.generBasisRed(objU.MlBas.levelBas{1});
%              Ny = objV.MlBas.levelBas{1}.generBasisRed(objV.MlBas.levelBas{1});
%              plot(objU.MlBas.levelBas{1}.plotVector,NU)

             % Find control points.
             U = objU.MlBas.levelBas{1}.knotVector;
             p = objU.MlBas.levelBas{1}.p;
             a = objU.MlBas.levelBas{1}.a;
             b = objU.MlBas.levelBas{1}.b;
             npts = length(U)-p-1;
             tauTest=linspace(a,b,npts);
             nstep=p;
             colmatTest=spcol(U,p+1,brk2knt(tauTest,nstep)); 
             [i1,i2]=size(colmatTest);      
             Nx = colmatTest(1:nstep:i1,:);
             DNx = colmatTest(2:nstep:i1,:);
             ctrlpts = Nx \ tauTest';
             ic=0;
             for j=1:npts
                 for i=1:npts
                     ic = ic+1;
                     xl(1,ic) = ctrlpts(i);
                     xl(2,ic) = ctrlpts(j);
                 end
             end
             xl;
%              hold all
%              plot(xl(1,:),xl(2,:),'k*', 'markers',4)
%              hold all

             Stiffn = zeros(nOF); % Initialize stiffness matrix.
             rhs = zeros(nOF,1); % Initialize right hand side.
             for elv =  1: nOEV
                 elv;
                 for elu =  1: nOEU
                     elu;
                     % Evaluation of ACTIVE basis functions and its derivatives at element elu.
                     % sU contains all the Gauss points and wU contains their weights.
                     % lIndex cointains the level of each active function and their indices.
                     [bValU,gradValU,lIndexU,r,wr] = evalEl(objU.MlBas,elu);
                     [bValV,gradValV,lIndexV,s,ws] = evalEl(objV.MlBas,elv);
                     gauss = [r,s];
%                      hold all
%                      plot(r(:),s(:),'k*', 'markers',4)
                     
                     % Find evaluation of basis functions and its derivatives at gauss points.
                     nstep=p;
                     for i=1:length(r)
                        Test=spcol(U,p+1,brk2knt(r(i),nstep)); 
                        [size1,size2]=size(Test);      
                        Gx(i,:) = Test(1:nstep:size1,:);
                        Gy = Gx;
                        DGx(i,:) = Test(2:nstep:size1,:);
                        DGy = DGx;
                     end
                     
                     % Find the tensor product of them.
                     for k=1:length(r)
                         ic=0;
                         for j=1:size(Gx,2)
                             for i=1:size(Gx,2)
                                 ic=ic+1;
                                 GxGy(k,ic)= Gx(k,i)*Gy(k,j);
                                 DGxGy(k,ic)= DGx(k,i)*Gy(k,j); 
                                 GxDGy(k,ic)= Gx(k,i)*DGy(k,j);
                             end
                         end
                     end
                     
                     % Find the Jacobian matrix dx/dξ,dy/dη.
                     for q=1:length(r)
                         for i=1:2
                             jacob(i,1,q)=0;
                             for k=1:size(GxGy,2) % Start loop for every basis function (and derivative).
                                 jacob(i,1,q)=jacob(i,1,q)+xl(i,k)*DGxGy(q,k);
                             end
                             jacob(i,2,q)=0;
                             for t=1:size(GxGy,2)
                                 jacob(i,2,q)=jacob(i,2,q)+xl(i,t)*GxDGy(q,t);
                             end
                         end
                     end
                     jacob;
                     
                     % Find the determinant of Jacobian matrix and construct its inverse matrix.
                     for k=1:length(r)
                         detjacob(k)=jacob(1,1,k)*jacob(2,2,k)-jacob(1,2,k)*jacob(2,1,k);
                         invjacob(1,1,k)=jacob(2,2,k)/detjacob(k);
                         invjacob(2,2,k)=jacob(1,1,k)/detjacob(k);
                         invjacob(1,2,k)=-jacob(1,2,k)/detjacob(k);
                         invjacob(2,1,k)=-jacob(2,1,k)/detjacob(k);
                     end
                     detjacob;         
                     % Find the local Jacobian matrix dξ/dr,dη/ds.
                     for k=1:length(r)
                         detjlocal(k) = 1/4*detjacob(k);
                     end
                     detjlocal;
                     
                     % Find the tensor product of evaluation of active basis functions at gauss points.
                     for k=1:length(r)
                         ic=0;
                         for j=1:size(bValV,2)
                             for i=1:size(bValU,2)
                                 ic=ic+1;
                                 % Nu*Nv
                                 shploc3(k,ic)= bValU(k,i)*bValV(k,j);
                                 % DuNu*Nv
                                 shploc1(k,ic)= gradValU(k,i)*bValV(k,j); 
                                 % Nu*DvNv
                                 shploc2(k,ic)= bValU(k,i)*gradValV(k,j);
                             end
                         end
                     end

                     % Store the level and the index of each tensor basis function.
                     j=size(lIndexU,2);
                     for i = 1:j
                         lIndex(1,i) = lIndexU(1,i);
                         lIndex(2,i) = lIndexU(2,i);
                         for k = 1:j-1
                             lIndex(1,i+k*j) = lIndexU(1,i);
                             lIndex(2,i+k*j) = lIndex(2,i);
                         end
                     end

                     j=size(lIndexV,2);
                     for i = 1:j
                         lIndex(3,i) = lIndexV(1,i);
                         lIndex(4,i) = lIndexV(2,1);
                         for k = 1:j-1
                            lIndex(3,i+k*j) = lIndexV(1,i);
                            lIndex(4,i+k*j) = lIndex(2,k+1);
                         end
                     end                     
                     lIndex;
                     
                     % Find the global derivative.
                     for j=1:length(r)
                         for i=1:size(shploc3,2)
                             prod3(j,i) = shploc3(j,i);
                             tp(j) = shploc1(j,i)*invjacob(1,1,j)+shploc2(j,i)*invjacob(2,1,j);
                             prod2(j,i) = shploc1(j,i)*invjacob(1,2,j)+shploc2(j,i)*invjacob(2,2,j);
	                         prod1(j,i)=tp(j);
                         end
                     end
                                          
                     % Assemble local Stiffness matrix.
                     elRhs = zeros(size(prod3,2),1); % Initialize local right hand side.
                     elRInd = cell(size(prod3,2),1);
                     elSInd = cell(size(prod3,2));
                     elStiff = zeros(size(prod3,2)); % Initialize local stiffness matrix.
                     for ii0 = 1 : size(prod3,2) 
                         for j = 1:length(s)
                             for i =1:length(r)
                                 elRhs(ii0) = elRhs(ii0)+wr(i)*ws(j)*detjlocal(j)*detjacob(j)*objU.f(r(i),s(j))*prod3(i,ii0);
                             end
                         end
                         elRInd{ii0} = lIndex(:,ii0);
						 for j = 1:length(s)
                             for i =1:length(r)
                                 elStiff(ii0,ii0) = elStiff(ii0,ii0)+wr(i)*ws(j)*detjlocal(j)*detjacob(j)*(prod1(i,ii0).^2+prod2(i,ii0).^2);
                             end
                         end
                         elSInd{ii0,ii0} = [lIndex(:,ii0) lIndex(:,ii0)];
                         for ii = ii0+1 : size(prod3,2)
                             for j = 1:length(s)
                                 for i =1:length(r)
                                     elStiff(ii0,ii) = elStiff(ii0,ii)+wr(i)*ws(j)*detjlocal(j)*detjacob(j)*(prod1(i,ii0).*prod1(i,ii)+prod2(i,ii0).*prod2(i,ii));
                                 end
                             end
                             elSInd{ii0,ii} = [lIndex(:,ii0) lIndex(:,ii)];
                             elStiff(ii,ii0) = elStiff(ii0,ii);
                             elSInd{ii,ii0} = elSInd{ii0,ii};
                         end
                     end
                     elRhs;
                     elStiff;
                                          
                     % Assemble global Stiffness matrix.
                     
                     for l = 1 : objU.MlBas.level
                         ic = 0;
                         for kv = 1:length(objV.MlBas.levelBas{l}.activeIndex)
                             for k = 1:length(objU.MlBas.levelBas{l}.activeIndex)
                                 ic = ic+1;
                                if(l>1)
                                   ind_1 = 0;  
                                   test = 0;
                                   for tInd = 1 : l-1
                                       ind_1 = ind_1 + length(objU.MlBas.levelBas{tInd}.activeIndex)^2;
                                       test = test + length(objU.MlBas.levelBas{tInd}.activeIndex);
                                   end                                   
                                   ind_1 = ind_1 + ic;
                                   test = test + k;
                                else
                                   ind_1 = ic;
                                end
                             
                             for ii1 = 1 : size(prod3,2)                                 
%                                  if(elRInd{ii1} == [l;objU.MlBas.levelBas{l}.activeIndex(k);l;objV.MlBas.levelBas{l}.activeIndex(kv)])
                                 if((elRInd{ii1}(1) == l)&(elRInd{ii1}(2) ==objU.MlBas.levelBas{l}.activeIndex(k)))
                                     rhs(ind_1) = rhs(ind_1) + elRhs(ii1);
                                     if l == 1
                                         iLvl(k) = l;
                                         iBasisFctIndU(k) = objU.MlBas.levelBas{l}.activeIndex(k);
                                         iBasisFctIndV(k) = objV.MlBas.levelBas{l}.activeIndex(k);
                                     else
                                         iLvl(test) = l;
                                         iBasisFctIndU(test) = objU.MlBas.levelBas{l}.activeIndex(k);
                                         iBasisFctIndV(test) = objV.MlBas.levelBas{l}.activeIndex(k);
                                     end
                                 end
                             end
                             
                             for ll = l : objU.MlBas.level
                                 ic2 = 0;
                                 for kkv = 1:length(objV.MlBas.levelBas{ll}.activeIndex)
                                     for kk = 1:length(objU.MlBas.levelBas{ll}.activeIndex)
                                         ic2 = ic2+1;
                                        if(ll>1)
                                           ind_2 = 0;
                                           for tInd = 1 : ll-1
                                               ind_2 = ind_2 + length(objU.MlBas.levelBas{tInd}.activeIndex)^2;
                                           end
                                           ind_2 = ind_2 + ic2;
                                        else
                                           ind_2 = ic2;
                                        end

                                     for ii = 1 : size(prod3,2)
                                         for iii = ii : size(prod3,2)
%                                              if(elSInd{ii,iii} == [l ll; objU.MlBas.levelBas{l}.activeIndex(k)...
%                                              objU.MlBas.levelBas{ll}.activeIndex(kk);l ll;objV.MlBas.levelBas{l}.activeIndex(kv)...
%                                              objV.MlBas.levelBas{ll}.activeIndex(kkv)])
                                             
                                             if((elSInd{ii,iii}(1,1) == l)&(elSInd{ii,iii}(1,2)== ll)&...
                                                (elSInd{ii,iii}(2,1) == objU.MlBas.levelBas{l}.activeIndex(k))&...
                                                 (elSInd{ii,iii}(2,2)== objU.MlBas.levelBas{ll}.activeIndex(kk)))
                                                 Stiffn(ind_1,ind_2) = Stiffn(ind_1,ind_2) + elStiff(ii,iii);
                                                 Stiffn(ind_2,ind_1) = Stiffn(ind_1,ind_2);
                                             end
                                         end
                                     end
                                     end
                                 end
                             end
                             end
                         end
                     end
                 end
             end
             Stiffn;
             Stiffn(:,all(Stiffn==0))=[];
             Stiffn(~any(Stiffn,2),:)=[];
             rhs;
             rhs(~any(rhs,2),:)=[];
             objU.iLvl = iLvl;
             objU.iBasisFctInd = iBasisFctIndU;
             objV.iBasisFctInd = iBasisFctIndV;
         end
         
         function y = solveSyst(objU,objV,Stiffn,rhs,dBC)
             [sdLvl, sdInd] = objU.MlBas.getActiveFctIndU(objU.MlBas.levelBas{1}.a);
             sEval = objU.MlBas.evalBasis(objU.MlBas.levelBas{1}.a);
             sIndEval = find(sEval);
             sFLvl = find( sdLvl(sIndEval(1)) == objU.iLvl);
             sFInd = find(sdInd(sIndEval(1)) == objU.iBasisFctInd);
             sInd = intersect(sFLvl,sFInd);
             
             [edLvl, edInd] = objU.MlBas.getActiveFctIndU(objU.MlBas.levelBas{1}.b-objU.MlBas.levelBas{1}.resol);
             eEval = objU.MlBas.evalBasis(objU.MlBas.levelBas{1}.b-objU.MlBas.levelBas{1}.resol);
             eIndEval = find(eEval);
             eFLvl = find( edLvl(eIndEval(end)) == objU.iLvl);
             eFInd = find(edInd(eIndEval(end)) == objU.iBasisFctInd);
             eInd = intersect(eFLvl,eFInd);
             
             if(dBC.mixed == true)
                nOFU = objU.MlBas.nOF; % Number of basis functions in U direction.
                nOFV = objV.MlBas.nOF; % Number of basis functions in V direction.
                nOF = nOFU*nOFV; % Number of basis function in the tensor product U*V.
                A = zeros(nOF +1);
                b = zeros(nOF +1,1);
                if(strcmp(dBC.east,'Dirichlet') && strcmp(dBC.west,'Neumann'))
                    A(1,eInd+1) = 1;
                    A(eInd+1,1) = 1;
                    b(1,1) = dBC.eVal;
                    A(2:end,2:end) = Stiffn;
                    b(2:end) = rhs;
                    u = A\b;
                    y = u(2:end);
                elseif(strcmp(dBC.west,'Dirichlet') && strcmp(dBC.east,'Neumann'))
                    A(1,sInd+1) = 1;
                    A(sInd+1,1) = 1;
                    b(1,1) = dBC.wVal;
                    A(2:end,2:end) = Stiffn;
                    b(2:end) = rhs;
                    u = A\b; 
                    y = u(2:end);
                end
             elseif(strcmp(dBC.east,'Dirichlet'))
                nOFU = objU.MlBas.nOF; % Number of basis functions in U direction.
                nOFV = objV.MlBas.nOF; % Number of basis functions in V direction.
                nOF = nOFU*nOFV; % Number of basis function in the tensor product U*V.
                
%                 Set boundary conditions on these nodes.
%                 X X X X X
%                 X - - - X
%                 X - - - X
%                 X - - - X
%                 X X X X X
%                 
%                 for i = 1:eInd
%                     bounds(i) = i;
%                     if i == sInd
%                        for j = 1:eInd-1
%                            bounds(i+j*eInd-((eInd-2)*j-(eInd-2))) = i+j*eInd;
%                            bounds(i+j*eInd-((eInd-2)*j-(eInd-2))+1) = i+j*eInd+eInd-1;
%                        end    
%                     end
%                     for k = 0:eInd-1
%                         bounds(3*(eInd-1)+k) = 1+(eInd-1)*eInd+k;
%                     end
%                 end      
%                 bounds;
                
                bounds=[];
                bounds(1) = 1;
                for i = 1:eInd
                    bounds(2*i) = i*eInd;
                    if i < eInd
                       bounds(2*i+1) = bounds(2*i)+1;
                    end
                end  
                bounds;

                A = zeros(size(Stiffn,1) +length(bounds));
                b = zeros(size(Stiffn,1) +length(bounds),1);
                
                for i = 1:length(bounds)
                    A(i,bounds(i)+length(bounds)) = 1;
                    A(bounds(i)+length(bounds),i) = 1;
                end
                 
                for i = 1:length(bounds)
                    b(i,1) = dBC.wVal;
                end
 
                A(length(bounds)+1:end,length(bounds)+1:end) = Stiffn;
                b(length(bounds)+1:end) = rhs;
                A;
                b;
                u = A\b;
                y = u(length(bounds)+1:end);   

             elseif(strcmp(dBC.west,'Neumann'))
                y = Stiffn\rhs;
             end
         end
         
         function uhU = generSolThb(objU,objV,y)
            DU = [];
            DV = [];
            for ll = 1 : objU.MlBas.level
                CU = zeros(objU.MlBas.levelBas{ll}.sP,length(objU.MlBas.levelBas{ll}.activeIndex));
                CV = zeros(objV.MlBas.levelBas{ll}.sP,length(objV.MlBas.levelBas{ll}.activeIndex));
                for l = 1:length(objU.MlBas.levelBas{ll}.activeIndex)
                    basisFU = thbSplBasFun(objU.MlBas.levelBas{ll}.activeIndex(l),objU.MlBas,ll);
                    basisFV = thbSplBasFun(objV.MlBas.levelBas{ll}.activeIndex(l),objV.MlBas,ll);
                    CU(:,l) = basisFU.generOneBasisFun;
                    CV(:,l) = basisFV.generOneBasisFun;
%                     Z = CU(:,1)*CV(:,1)';
                end
%                 DU = [DU Z(1:2^(ll-1):size(Z,1),:)];
%                 DU = [DU CU(1:2^(ll-1):size(CU,1),:)];
%                 DV = [DV CV(1:2^(ll-1):size(CV,1),:)];
            DU = [DU CU(1:2^(ll-1):size(CU,1),:)];
            end
            uhU = y(1)*DU(:,1);
            for k = 2 : objU.MlBas.nOF
                uhU = uhU + y(k)*DU(:,k);
            end
%             ic = 0;
%             for j=1:size(DV,2)
%                 for i=1:size(DU,2)
%                     ic=ic+1;
%                     DUDV(:,ic)= DU(:,i).*DV(:,j);
%                 end
%             end
%             nOF = objU.MlBas.nOF*objV.MlBas.nOF;
%             uhU = y(1)*DUDV(:,1);
%             for k = 2 : size(y,1)
%                 uhU = uhU + y(k)*DUDV(:,k);
%             end
         end
    
      end
end
                                                                            