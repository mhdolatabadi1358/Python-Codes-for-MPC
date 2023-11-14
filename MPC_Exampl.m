%We run the code presented in:

%  https://it.mathworks.com/help/mpc/ug/solve-custom-mpc-quadratic-programming-problem-and-generate-code.html

%  explained in page 12 of the "Lecture notes" by Professor Mark Cannon's lecture notes for the Model Predictive 
%       Control class at University of Oxford

%  https://markcannon.github.io/assets/downloads/teaching/C21_Model_Predictive_Control/C21_MPC_Lecture_Notes.pdf 

clear all

A = [1.1 2; 0 0.95];
B = [0; 0.0787];
C = [-1 1];
D = 0;
Ts = 1;
sys = ss(A,B,C,D,Ts);
x0 = [0.5;-0.5]; % initial states at [0.5 -0.5]
Qy = 1;
R = 0.01;
Qb=C'*Qy*C; Nb=C'*Qy*D; Rb=D'*Qy*D+R;
K = -lqr(sys,Qb,Rb,Nb);

t=0:40;
lt=length(t);

x=x0;
Y_lqr=zeros(lt,1);
U_lqr=zeros(lt,1);
for i=1:lt
    Y_lqr(i)=C*x;
    u= K*x;
    u= max(-1,u); u=min(1,u);
    U_lqr(i,1)=u;
    x = A*x+B*u;        
end


Q = C'*C;
AA=(A+B*K)';
QQ=Q+K'*R*K;

Q_bar = dlyap(AA, QQ);

%----------------------------- formula (2.2) page 16
n=size(A,1);
m=size(B,2);
l=size(C,1);
N=16;  %---- Horizon
M=zeros(N*n,n);
for i=1:N; M((i-1)*n+1:i*n,:)=A^i; end

CONV=zeros(N*n,N*m);
for i=1:N; for j=1:i; CONV((i-1)*n+1:i*n,(j-1)*m+1:j*m)=(A^(i-j))*B; end;end

%----------------------------- formula at page 17
Q_hat=zeros(n*N,n*N);
for i=1:N-1; Q_hat((i-1)*n+1:i*n,(i-1)*n+1:i*n)=Q; end
Q_hat((N-1)*n+1:N*n,(N-1)*n+1:N*n)=Q_bar;
R_hat=zeros(l*N,l*N);
for i=1:N; R_hat((i-1)*l+1:i*l,(i-1)*l+1:i*l)=R; end
  
H = CONV'*Q_hat*CONV + R_hat;  %---- formula (2.4) page 17
F = CONV'*Q_hat*M;

H=(H+H')/2; %---------- to make the matrix symmetric

Ac = [eye(N);-eye(N)];  b0 = ones(2*N,1);  %---- formula (2.12) page 27

x = x0; 
Y_MPC=zeros(lt,1);U_MPC=zeros(lt,1);

counter=1;
for ct = 1:lt
    Y_MPC(ct,1)=C*x;
    u1 = quadprog(H,F*x,Ac,b0);   %--------------  By Mohan
    u = u1(1:m,1); 
    U_MPC(ct,1)=u;
    x = A*x+B*u;    
end

figure
subplot(2,1,1)
stairs(t,U_MPC(1:lt),'linewidth',1.2)
hold on
stairs(t,U_lqr(1:lt),'--','linewidth',1.2);
xlabel('time')
ylabel('u')
ylim([-1.5 3])
subplot(2,1,2)
plot(t,Y_MPC(1:lt),'linewidth',1)
hold on
plot(t, Y_lqr(1:lt),'--','linewidth',1)
xlabel('time')
ylabel('y')
ylim([-4 6])
legend('Constrained MPC')

