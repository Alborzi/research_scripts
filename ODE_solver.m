%***************************MAIN MATLAB Script****************************
%--------------------------------------------------------------------------
%Reading the thermodynamics file
[no a b c H S  spec]=textread('AFTS.txt','%f %f %f %f %f %f %s',-1);
%Reading the chemical mechanism file
[ar br cr dr er fr A_a n_exp E_a]=textread('mechanism.txt','%f %f %f %f %f %f %f %f %f',-1);
%--------------------------------------------------------------------------
%input
T=423;                                %     T of the Reaction
RGAS=1.987e-3;                        %     Uinversal Constant of Gases
n=length(ar);                         %     No of reactions  
m=19;                                 %     No of species
%--------------------------------------------------------------------------
%construction of Stoichiometric matrix based on reaction mechanism
alpha=zeros(n,m);
for ii=1:n
        
        if br(ii)~=0 
           alpha(ii,abs(br(ii)))=alpha(ii,abs(br(ii)))+abs(br(ii))/br(ii);
        end
        if cr(ii)~=0
            alpha(ii,abs(cr(ii)))=alpha(ii,abs(cr(ii)))+abs(cr(ii))/cr(ii);
        end
        if dr(ii)~=0 
            alpha(ii,abs(dr(ii)))=alpha(ii,abs(dr(ii)))+abs(dr(ii))/dr(ii);
        end
        if er(ii)~=0
            alpha(ii,abs(er(ii)))=alpha(ii,abs(er(ii)))+abs(er(ii))/er(ii);
        end
        if fr(ii)~=0
            alpha(ii,abs(fr(ii)))=alpha(ii,abs(fr(ii)))+abs(fr(ii))/fr(ii);
        end
               
end

%--------------------------------------------------------------------------


       for ii=1:length(ar)
           %k_forward(ii)=A_a(ii)*exp(-E_a(ii)/(RGAS*T));
           k_forward(ii)=A_a(ii)*(T.^n_exp(ii))*exp(-E_a(ii)/(RGAS*T));
                  end 
          
                ks=[k_forward' ]; 
                
 %------------------------------------------------------------------------
 %---------------------- Experimental data -----------
 experimental=[

 ];
 
   texp=experimental(:,1)*60;
   Cexp=experimental(:,2);
 %-------------------------------------------------------------------------
   
 %----------------------Specify chemical species at time=0---------------
 % The order of species is the same as the species index in the AFTS.txt
   
C0=zeros(19,1); % the first number indictes the total number of chemical species in the mechanism
C0(2)=;        
C0(4)=;         
C0(5)=;              
%-----------------------------------------------------------------------------
%---------------------------------------------------------------------------
  t=0:0.1:1500;  %Time for numerical integration
  
options=odeset('AbsTol', 1e-20,'RelTol',1.e-10);
[t C]=ode15s(@dCdt,t,C0,options,ks,alpha,br,cr);
CAFTS=C;
timeAFTS=t;
hold all
plot(t/60,C(:,2),texp/60,Cexp,'o');
%plot(t,C(:,2),texp2,Cexp2,'x');
%************************************************************************************************************
%******************************************* Function for reaction rate*************************************

% The input and output arguments are in the order required for use with the
% Matlab main ODE solvers.
% Inputs:
% t = time
% Cs = concentrations at time t
% ks = rate coefficients
% alpha = matrix of stoichiometric coefficients
% Cs and ks are both COLUMN vectors
% Output:
% Cp = derivative of concentration with respect to time
%
[nrxn nspec]=size(alpha);

%%%%%%%%%%%%%%%%%%%%%%%%
%******Function*********
%%%%%%%%%%%%%%%%%%%%%%%%
for i =1:nrxn;
    if cr(i)==0
        rf(i)=rf(i)*C(abs(br(i)))^(-alpha(i,abs(br(i))));
    elseif br(i)==0
        rf(i)=rf(i)*C(abs(cr(i)))^(-alpha(i,abs(cr(i))));
    else
        rf(i)=ks(i,1)*C(abs(br(i)))*C(abs(cr(i)));
    end
end
r=rf;
Cp=alpha'*r;

return
%*****************************************************************************************************************
***************************************** Species list""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
% Columns 2-6 are place-holders for coeficients for polynomials(thermochemistry and Kequlibrium calculation) 
% However, at the moment they are not read by code. The user only need to simply specifiy the name of the species in 
% the last column. In case of any problems please contact e.alborzi@sheffield.ac.uk
%-------------------------------------------------------------------------------------------------------------
1       3.9000          0.23680         -0.000084           -69.66          155.02      R.
2   	6.8760       	0.00044         0.0000031           0.000           69.074      O2		
3       21.794         	0.03262     	0.0000000           21.019	        106.99      RO2.
4       0.1500          0.25010         -0.000100           -64.53          134.91      RH
5   	7.1500          0.23690         -0.000100           -71.00          169.87      RO2H
6   	7.6950          0.00250         0.0000000           9.291           149.99      R2
7   	66.180          4.19000         8267.7000           -116.53         215.33      RO.	
8   	9.9700          0.23240         -0.000100           -91.06          169.74      OH.
9   	7.9710          -0.0009         0.0000000           -145.16         288.62      ROH
10  	16.190      	0.23110         -0.000100           -57.76          45.070      Rprime
11  	7.7000          0.50002         -0.000200           215.53          -108.67     Carbonyl
12      7.7000          0.50002         -0.000200           215.53          -108.67     H2O
13      7.7000          0.50002         -0.000200           215.53          -108.67     alkane
14 		7.7000          0.50002         -0.000200           215.53          -108.67		ROOOOR
15 		7.7000          0.50002         -0.000200           215.53          -108.67		ROR
16   	7.7000          0.50002         -0.000200           215.53           -108.67    H.
17 		7.7000          0.50002         -0.000200           215.53           -108.67    HOO.
18      7.7000          0.50002         -0.000200           215.53           -108.67 	ROOR
19      7.7000          0.50002         -0.000200           215.53           -108.67 	Carbonyl2
%*****************************************************************************************************************
%***************************************** Chemical kinetic mechanism********************************************
% By convention reactant species have the negative sign and products have positive sign. The species index are taken from the species 
% list.
%-------------------------------------------------------------------------------------------------------------

1		 0		-4	    16      1       0       5.4E11      0       61.5
2       -4		-2	    1      17       0 		2.7E7    	0       55.1 
3       -1      -2      3   	0   	0       7.5E8       0    	0.000
4   	-3      -4  	5       1   	0       4.75E10     0   	14.45			
5   	-1      -1  	6   	0       0   	2.58E11     0	    0.000	
6       -5      -0  	7   	8   	0       3E13        0    	33.0
7		-5      -5      7       3       12      1E10        0   	24.1
8       -3      -7      11       5      0        1E9        0       10.60 
9       -8      -4      12  	1       0       2.5E9       0	    5
10      -7       0      10      11      0       2.8E12      0	    9.0		
11      -7      -4      9       1       0       1.85E9      0   	5.3	
12      -10  	-4      13  	1   	0   	1E9         0   	10.10	
13      -3       0   	1       2  	    0   	1E16        0   	19.10
14      -7      -7		11		9		0		1E9	        0		17.5
15       -3		-3      14      0       0       9.6E10		0		0
16      -14     0		7		7		2		2.6E13      0       0.8 
17       -7     -1      11		4		0		2E9 		0		7.5
18       -7     -1      19      4       0       3.4E9       0       9.8
19       -7    -1      15       0       0       2.7E9       0       0
20       -7    -7      18       0       0       2.7E9       0       0    





