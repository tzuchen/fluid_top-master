% The is the matlab routine to test the send and receive 
% function between petsc and matlab
%
%
%
%
%   Tzu-Chen Liang

addpath /home/tzuchen/snopt7111/snopt-dev/matlab
addpath /home/tzuchen/petsc/petsc-2.3.0/bin/matlab
addpath /home/tzuchen/proj/fluid_top/trunk/matlab_send
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 16;
n = 5;
x_send = rand(1,1);
v_send = rand(m*n,1);  

p_send = speye(50);
B = diag(ones(m-1,1),2)+diag(ones(m+1,1));;
B = B(:,1:end-1);
B(2,1)=1;
B = rand(1000,1000);
p_send = sparse(B);

PetscInitialize(5,'testcommmatlab');
%a = openport(5050);
a = openport();

disp(sprintf('\n\nExample :Send and receive a scalar'))

    send(a,x_send);
       disp('x has been sent')
    x_receive = receive(a);
       disp('x has been received!')
    x_receive = x_receive(1);  % ScalarView send "size" copies of x back! 

       disp(sprintf('The difference between x_send and x_receive is %d',norm(x_send-x_receive)))

disp(sprintf('\n\nExample :send and receive a SEQ vector!'))
  
    send(a,v_send);
       disp('v has been sent!')
    v_receive = receive(a);
       disp('v has been received!')
    v_receive = v_receive(1:m*n);
 
       disp(sprintf('The difference between v_send and v_receive is %d',norm(v_send-v_receive)))

disp(sprintf('\n\nExample :send and receive a MPI vector!'))
    
    send(a,v_send);
       disp('v has been sent!')
    v_receive = receive(a);
       disp('v has been received!')

  disp(sprintf('The difference between v_send and v_receive is %d',norm(v_send-v_receive)))


disp(sprintf('\n\nExample :send an MPI Sparse matrix!'))
   
    send(a,p_send);
        disp('p has been sent!')
    %p_receive = receive(a);
    %    disp('p has been received!')





pause(5)




%send(a,p)

closeport(a)
