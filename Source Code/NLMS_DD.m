%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Digital Signal Processing                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          As a Part of Digital Signal Processing Laboratory              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Acoustic Echo Cancellation with Double Talk Detection            %
%       Done By: 1) Suraj Khan (22455660)                                 %
%                2) Shrishti Saha Shetu (22464279)                        %
%                3) Raza Azam ()                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input values:
%   x = loudspeaker signal
%   d = microphone signal
%   N = Filter Length
%   alpha = step size
%   mem_time = Number of samples to take while calculating the PSD
%   freeze_time = Seconds for which the filter to turn off after detection
%   del_processing_time = Seconds by which the filtering is delayed
%
%return values:
%   e = estimation error (error signal) for each iteration
%   h_hat = estimate of h after the last iteration

function [e, h_hat, d_hat,h_str] = NLMS_DD(x, d, N, alpha, fs, mem_time, freeze_time, del_processing_time,h_str)

    %regularization constant to avoid divison by zero
    DELTA = 0.000001;

    %length of signals
    Lx = length(x);

    %init variables
    xk = zeros(N,1);        %last N samples of x
    
    N2=mem_time*fs;
    xk1 = zeros(N2,1);        %last N samples of d
    dk1 = zeros(N2,1);        %last N samples of d
    h_hat = zeros(N,1);
    e = zeros(Lx,1);
    d_hat =zeros(Lx,1);
    p=0;
    
    end_index= -freeze_time*fs;

    
   for k = 0 : Lx-1
        xk1 = [x(1+k); xk1(1:N2-1)];
        dk1 = [d(1+k); dk1(1:N2-1)];
         
       if(~(mod(k,1000)))       % Compute PSD after every 1000 Samples

            %Compute PSD
            cross_psd=cpsd(xk1,dk1);
            acf_mic=pwelch(dk1);
            acf_ref=pwelch(xk1);
            
            msc=((cross_psd')*cross_psd)./((acf_ref.')*(acf_mic));

            %Check Activity
            if(msc<0.15)
                  end_index= k-del_processing_time*fs;
                  fprintf('\n Filter blocked at %f secs with new K at %f',k/fs,end_index/fs);
            end
       end
        
        if(k>=del_processing_time*fs)
             %Store the Last N samples
            xk = [x(1+p); xk(1:N-1)];
            % Compute the output 
            y_est = h_hat' * xk;
            %Compute the Error
            e(p+1) = d(p+1) - y_est;
             %Update Coffiencients only when the sample is less than the freeze
             %index
            if (~(p<=((end_index-del_processing_time*fs)+freeze_time*fs)))
                  h_hat = h_hat + alpha *( e(p+1).* xk / (xk'*xk + DELTA));
            end
            p=p+1; 
        end
       % for i=1:2001
            h_str(k+1,:)=[h_hat.'];
      %  end
   end       
end
