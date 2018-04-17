%% This function is plotting the the share of the different reflection and absorption processes
% Only the first of the simulated batches will be displayed

fig1=figure
area(lambda,[A_si_1, R_ext_1, R_esc_1, A_grid_1, A_back_tot_1, A_AR_1, A_FCA_1])
axis([300 1200 0 1])

L=legend('A_{si}', 'R_{ext}', 'R_{esc}', 'A_{grid}', 'A_{back}', 'A_{AR}', 'A_{FCA}','location','South')
L.TextColor='w';
legend boxoff