function [adja1,adja2] = prep_data_for_GLEE(P1, Q1, Eg1 ,P2, Q2, Eg2)
    adja1 = genAdjaMatrix(P1, Q1, Eg1);
    adja2 = genAdjaMatrix(P2, Q2, Eg2);

end

