function [A, B, C, D] = ConvDiscSISO(Ac,Bc,Cc,Dc,Ts,IN,OUT)

%converts any system to a discrete time SISO system.  The converted system
%does NOT nessesarily have the same behavior as the original system, only
%use for testing purposes

%Input D is assumed to be 0

%Inputs are any stable system matricies
%OPTIONAL INPUTS:
%    IN: input dimention
%    OUT: Output dimension
%Outputs are a related system that is discrete time and SISO



Ac = full(Ac);
Bc = full(Bc);
Cc = full(Cc);
%check if we need to convert to SISO or to a discrete time system
[~, num_in] = size(Bc);
[num_out, ~] = size(Cc);
CTsys = max(abs(eigs(Ac))) > 1;

needToChange = CTsys || num_in > 1 || num_out > 1;

%if 6 inputs are supplied, we select the requested input-output dimentions.
%IF 4 inputs are supplied, default to input_dim = output_dim = 1
if nargin == 6
    input_dim = IN;
    output_dim = OUT;
else
    input_dim = 1;
    output_dim = 1;
end


if needToChange
    %convert to discrete time SISO system.
    if num_in > 1
        Bc = Bc(:,input_dim);
    end
    if num_out > 1
        Cc = Cc(output_dim,:);
    end
    if num_out > 1 || num_in > 1
        fprintf('Changed to SISO (%d:%d)\n',input_dim,output_dim)
    end
    Dc = 0;
    if CTsys
        sysc = ss(full(Ac),Bc,Cc,Dc);
        sysd = c2d(sysc,Ts);
        A = sysd.A;
        B = sysd.B;
        C = sysd.C;
        D = sysd.D;
        fprintf('Converted from continuous to discrete\n')
    else
        A = Ac;
        B = Bc;
        C = Cc;
        D = Dc;
    end
else
    fprintf('Did not change the system\n')
    A = Ac;
    B = Bc;
    C = Cc;
    D = Dc;
end