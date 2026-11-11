function u_exact = exact_sol(x,y)
u1=sin(pi*x)*sin(pi*y);
u2=x*(x-1)*y*(y-1);
u_exact=[u1;u2];
end