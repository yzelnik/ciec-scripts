function dy = ode_glv(~,y,kv,am)

dy = y(:).*(kv + y(:)'*am')';

end