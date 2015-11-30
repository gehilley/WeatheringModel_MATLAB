function [LBC,RBC] = bt_boundaryConditions(u,uo)

LBC = repmat(reshape(uo,1,length(uo)),2,1);

RBC = repmat(u(end,:),2,1);
