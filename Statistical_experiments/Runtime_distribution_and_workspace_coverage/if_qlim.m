function f=if_qlim(q,qlim)
% If joint angle is in its motion range?
if all(q>=qlim(:,1) & q<=qlim(:,2))
    f=1;
else
    f=0;
end