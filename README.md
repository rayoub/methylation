This code is based on Martin Sill's code from here https://github.com/mwsill/mnp_training.
I just modularized the code and dispensed with validation parts they used for their paper. 
Also, I upgraded from randomForest package implementation of the Random Forest algorithm to the ranger package. 
This makes training much faster and uses multiple cores internally. 
