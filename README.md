# predictivePursuit
Code accompanying the paper:

Kashyap, H. J., Detorakis, G., Dutt, N., Krichmar, J. L., & Neftci, E. (2018). A Recurrent Neural Network Based Model of Predictive Smooth Pursuit Eye Movement in Primates. In International Joint Conference on Neural Networks (IJCNN).

http://www.socsci.uci.edu/~jkrichma/Kashyap-PredicitvePursuit-IJCNN2018.pdf

The code is tested using Matlab 2017a. If you use this code in your research, please cite the above paper.

Files:

main_sp_predict.m (main script to run the prediction and the initiation experiemnts)

sp_predict_ramp_pertur.m (script to run the predictive pursuit experiment with perturbation and phase shift)

sp_predict_random_target.m (script to run the predictive pursuit experiemnt for a random target)

get_accelaration.m (Calculates pursuit accelaration during initiation)

get_RS.m (plots retinal slip)

getWeights.m (RNN synapse weights)
