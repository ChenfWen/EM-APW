# EM-APW
Reproduction code for Economic Modelling Paper：
"Antipollution willingness and urban population distribution: Evidence from Chinese cities"

The files include two categories: Stata files and MATLAB files


The original data is 'data_od_vc_fdi_ind.csv' in the 'Stata' folder. The following is the implication of indices in the original data:

  city ---> the name of the city
  
  areacode ---> the code of the city
  
  year ---> year
  
  pop ---> population of the city
  
  gdp ---> gross domestic production of the city
  
  pm2.5 ---> the concentration of PM2.5 of the city
  
  labor ---> the number of workers of the city
  
  capital ---> the capital of the city
  
  solid_cut ---> 
  
  solid ---> 
  
  pop_huji ---> the population registered in the city (born in this city)
  
  o_meanyear ---> the migration cost from this city to others
  
  d_meanyear ---> the migration cost from other cities to this city
  
  vc ---> the air speed over the city
  
  FDI ---> foreign direct investment
  
  areacode2 ---> the province code that the city belongs to
  
  bl ---> the coefficient of labor
  
  bc ---> the coefficient of capital
  
  bp ---> the coefficient of pollutant
  


Stata files is to calculate the two characteristics (Productivity and APW) and do the empirical regression. The following is the relationship between the code files and results used in the MATLAB code and manuscript:

  cal_result.do ---> data_ctrl.csv
  
  cal_result_vc.do ---> data_ctrl_vc.csv
  
  reg_result.do ---> Table 2 

MATLAB files is to impose the counterfactual experiments. The C*-like files are the functions, in which A, ame, tau and rho denote productivity, amenity, migration friction, and APW respectively. (For example, CA is the counterfactual function for eliminating the divergence of productivity among cities.) The following is the relationship between the code files and results in the manuscript:

  Robust_eps.m ---> Figure 2 and 8
  
  Counterfactual_one.m ---> Figure 3
  
  Counterfactual_three.m ---> Figure 4
  
  Externality_A.m ---> Figure 5
  
  Externality_ame.m ---> Figure 6
  
  Externality_nochar.m ---> Figure 7
  
  Counterfactual _education.m ---> Figure 9 and 10
  
  Robust_vc.m ---> Figure 12
  
  Robust_ind ---> Figure 11
