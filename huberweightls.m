function w=huberweightls(data, k);
% the weight function for huber's robust regression function.
%
%(c) Copyright Lingsong Zhang (lingsong@purdue.edu)

w=k./abs(data);
w(w>1)=1;


