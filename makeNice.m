function [plot2D] = makeNice(signal,shift)
%% makeNice(signal) rearranges the fourier transformed signal appropriately for the rest of the analysis code
    if shift 
        niceSignal = rot90(transpose(fftshift(signal)),2);
    else
        niceSignal = rot90(transpose(signal),2);
    end
plot2D.Absolute = abs(niceSignal);
plot2D.Real = real(niceSignal);
plot2D.Imaginary = imag(niceSignal);


end

