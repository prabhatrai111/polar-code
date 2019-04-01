%%% ECC PROJECT BY PRABHAT,KRISHNA %%%
%%% POLAR ENCODING %%%

clear all; close all; clc; 

N=128; K=64; Ec=1; N0=2;

% N  =  Blocklength
% K  =  Message length
% Ec =  Taking the BPSK symbol energy
% N0 =  Noise power spectral density (we are taking N0/2 = 1) 
% Polar code construction at the given SNR (taking default 0dB) 

n = log2(N); % For number of Stages

% Constructing using Bhattacharyya parameter
% SNR = (Ec/N0) = (SNR_OF_AWGN/2), given in dB
z = zeros( N , 1 );
SNRdB = 2 * Ec / N0;

% Actual Bhattacharya Parameter for AWGN CHANNEL = exp(-Ec/N0), we are
% taking in log domain
z(1) = -10 ^ ( SNRdB / 10 );

fprintf("tree for encoding by bhattacharya parameter \n");
for lev = 1 : n
             B = 2 ^ lev;
    for j = 1 : B/2
        Temp = z ( j );
        z( j ) = logdomain_diff( log( 2 ) + Temp, 2 * Temp );  %  2 z - z ^ 2
        z( B / 2 + j ) = 2 * Temp;                          %  z ^ 2
        fprintf("j=%d, B=%d, T=z(%d)=%f, z(%d)=%f,z(%d)=%f \n",j,B,j,Temp,j,z(j),B/2+j,z(B/2+j));
    end
end
fprintf("tree encoding by bhattacharya parameter is done \n");

% sorting in Ascending, 1 is to sort in columns
[~,indices] = sort( z, 1, 'ascend');  
lookup = zeros( N, 1 );
for j = 1 : K 
    lookup( indices( j ) ) = 1; % locations in z containing least K values
end

% BINARY Message input
input = rand( K, 1 ) > 0.5; 

% Mapping of inputs to 'K' ascending index
   % Encode 'K' message bits in 'input' and
   % Output will be 'N' bits with frozen included
   % lookup(i)==  0     == bit-i is a frozenbit
   % lookup(i)==  1      == bit-i is a messagebit

d = lookup;               % loads all bits into d
d( lookup == 1 ) = input; % -1's will get replaced by message bits

% Encoding with xor/modulo addition
fprintf("Polar encoding is starting \n");
for i = 1 : n
           B = 2 ^ ( n - i + 1);
         lev = 2 ^ ( i - 1 );
    for j = 1 : lev
        bse = ( j - 1 ) * B; 
        for l = 1 : B/2
            ppp = d( bse + l );
            % Modulo addition in every stage and level
            d( bse + l ) = mod( d( bse + l ) + d( bse + B/2 + l), 2 );
  fprintf("d(%d)=%d xor d(%d)=%d => d(%d)=%d\n",bse+l,ppp,bse+B/2+l,d(bse+B/2+l),bse+l,d(bse+l));
        end
    end
end
fprintf("Polar encoding is over \n");

% where d is the polar encoded output
d;

% Transmission of polar codes with BPSK modulation and AWGN Noise
y = ( 2 * d - 1 ) * sqrt( Ec );                  % BPSK modulation
bpsk_output = y + sqrt( N0 / 2 ) * randn( N, 1 ); % AWGN (addition) only real part 
 
 
 
