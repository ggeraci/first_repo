19 October 2013:

Bug found: the last row of matrix which has all 1's keep its 1's even for j+k!=0
Remedy 1: loop in reverse from bif j and k to small ones >>> AMG and BiCG work (we have to pin the last p, to keep the sparsity pattern)
Remedy 2: Define an array of matrices J and K >>>> The AMG does not work!, cannot have last row with all 1's! NOW: last row merely pin the last p value to zero.



20 October 2013:

Bug found: the RHS of RU and U was not wero for the last slice, so it changes the previous value for the characteristic boundary conidition we are solving.


20 October 2013:

Math Bug in poisson consistency equation:

int div.u dv = int u dA ===>> sum(div) * dx * dy * dz = (sum(u_o)-sum(u_i))* dy* dz


22 October 2013:

Bug found: I had RU_int=ru_new before updating ru_new with pressure gradient


27 October 2013:

Issues with particle and flow data saving/loading resolved. Now the restart files work. Verified for different cases.

3 Novemebr 2013:

Bug in SEND_RECV_BLOCKING_CUM_INOUT: the function it self was correct, but never called from update_ghost_cum!!