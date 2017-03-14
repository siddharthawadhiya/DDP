function residue = model (t,y,yp,rxobj,N,Er,G,I,J)

phi = y(1:N);
theta = y(N+1);
P = y(N+2);
Vr = y(N+3);
rho = y(N+4);

phip = yp(1:N);
thetap = yp(N+1);
Pp = yp(N+2);
Vrp = yp(N+3);
rhop = yp(N+4)

rs = rxobj(phi,theta,P);



resphi = phip - rs.ralpha;
resT = thetap - ((rs.ra' * rs.halpha) + (theta +1/rs.D) * sum(phip))/sum(phi);
resP = Pp + Er(1)*G*Vr + Er(2)*I*((Vr)^2);
resVr = J*(1/sum(phi))*(1/(theta +1/rs.D))*P*Vr*(1/t) + J*rhop*Vr + J*(1/sum(phi))*(1/(theta +1/rs.D))*P*Vrp ;
resrho = rho - P*(1/sum(phi))*(1/(theta +1/rs.D));
residue = [resphi; resT; resP; resV; resR];

end
        