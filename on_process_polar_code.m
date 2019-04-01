%%% ECC PROJECT BY PRABHAT,KRISHNA %%%
%%% POLAR ENCODING %%%

clear all; close all; clc; 
N=2048; K=N/2; Ec=1; N0=2;

% N  -  Blocklength % K  -  Message length (Rate = K/N)
% Ec -  Taking the BPSK symbol energy
% N0 -  Noise power spectral density (we are taking N0/2 = 1) 
% Polar code construction at the given design-SNR (default 0dB) 

n = log2(N); % For number of Stages

% Constructing using Bhattacharyya parameter
% SNR = (Ec/N0) = (SNR_OF_AWGN/2), given in dB
z = zeros( N , 1 );
SNRdB = 2 * Ec / N0;
SNR = 10^( SNRdB / 10 );

% Actual Bhattacharya Parameter for AWGN CHANNEL = exp(-Ec/N0)
z(1) = -SNR;
for lev = 1 : n
             B = 2 ^ lev;
    for j = 1 : B/2
        T = z ( j );
        z( j ) = logdomain_diff( log( 2 ) + T, 2 * T );  %  2 z - z ^ 2
        z( B / 2 + j ) = 2 * T;                          %  z ^ 2
    end
end

% sorting in Ascending, 1 is to sort in columns
[~,indices] = sort( z, 1, 'ascend');  
lookup = zeros( N, 1 );
for j = 1 : K 
    lookup( indices( j ) ) = - 1; % locations in z containing least K values
end

% BINARY Message input
ebn_db = 1 : 1 : 20;
ebn = 10.^(ebn_db/10);
for tt = 1 : length(ebn)
    e=0;   
    for ii = 1 : N
        input = rand( K, 1 ) > 0.5;  
        % Encode 'K' message bits in 'u' and
        % Return 'N' encoded bits in 'd'
        % lookup(i)==  0     ==> bit-i is a frozenbit
        % lookup(i)== -1     ==> bit-i is a messagebit
        d = lookup; % loads all bits into d, including -1
        d( lookup == -1 ) = input; % -1's will get replaced by message bits below
        d_look_input = d;
        for i = 1 : n
                   B = 2 ^ ( n - i + 1);
                 lev = 2 ^ ( i - 1 );
            for j = 1 : lev
                bse = ( j - 1 ) * B; 
                for l = 1 : B/2
                    % Modulo addition in every stage and level
                    d( bse + l ) = mod( d( bse + l ) + d( bse + B/2 + l), 2 );
                end
            end
        end
        % where d is the polar encoded input
        polar_encoded_d = d;
        % Transmission of polar codes with BPSK modulation and AWGN Noise
        x_bpsk = ( 2 * polar_encoded_d - 1 ) * sqrt( Ec );
        noise = sqrt( N0 / (2*ebn(tt)) ) * randn( N, 1 );
        received_bits = x_bpsk + noise; % BPSK + AWGN

        for i = 1 : N
            if received_bits(i) > 0.1
                d_hat(i) = 1;
            else
                d_hat(i) = 0;
            end
        end
         d_hardML_hat = d_hat; 

         for i = 1 : n
             B = 2 ^ ( n - i + 1);
             lev = 2 ^ ( i - 1 );
             for j = 1 : lev
                 bse = ( j - 1 ) * B; 
                for l = 1 : B/2
                    % Modulo addition in every stage and level
                    d_hardML_hat( bse + l ) = mod( d_hardML_hat( bse + l ) + d_hardML_hat( bse + B/2 + l), 2 );
                end
             end
         end                        
        d1 = lookup; % loads all bits into d, including -1 
        for i=1:N
            if d1(i) == -1
                polar_decoded(i) = d_hardML_hat(i);
            else
                polar_decoded(i) = 0;
            end
        end
        if polar_decoded == d_look_input.'
            e=e;
        else
            e=e+1;
        end
        ii = ii + 1;
    end
    er(tt)=e/N;
end
figure; semilogy(ebn_db,er);