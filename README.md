# Dissertation project - Utilising machine learning to investigate electrical activity in ventricular myocytes

Machine Learning is a passion of mine. The project used a novel machine-learning based approach to investigate the effect that changing multiple ion conductances has on the transition from normal cell electrical activity (action potentials) to abnormal cell electrical activity that can cause lethal arrhythmias (early afterdepolarizations).

For my project, while handing a heavily imbalanced dataset, I trained an excellent XGBoost classifier which had a misclassification rate of just 0.07%. I also used SHAP values (a concept from game theory) to rank the importance of ion conductances in making predictions.

The analysis also compared the performance of multiple machine learning algorithms (multivariable logisitic regression, decision trees, random forests and XGBOOST). and balancing techniques (undersampling, oversampling, SMOTE, ROSE). In the end, XGBoost was the winner with the added benefit of being compatible with the SHAPforxgboost package in R. 

Overall the study proved that machine learning is a great technique for exploring action potential dynamics. The XGBoost model is an excellent classifier for distinguishing between EADs and normal action potentials with Random Forests as a close second. Additionally, using AUC as a metric was enough to handle the class imbalance problem. 

The study also confirmed previous findings that an IKr block and sufficient inward ICa,L play important roles in EAD genesis and are in fact the most influential currents in EAD genesis. More generally, this study supports previous research that state that both a decrease in repolarizing currents and an increase in depolarizing currents promote EAD genesis


## Files
Dissertation data contains 8 matlab files. The endothelial O’Hara Rudy mathematical model was simulated by another MSc student to create the simulations in the matlab files. Altogether, there are 212,249 observations (simulations) and 11 variables. These variables are cycle length, APD90, class (the dependent variable) and 8 independent variables - GNa, Gto, PCa, GKr, GKs, GK1, GNCX, GpCa. These 8 independent variables are altered channel conductances/permeability used to predict the class variable (action potential class) and their values represent the scale factor at which they were altered.


## Dissertation can be found on my github page https://jhfran.github.io/mldissertation
