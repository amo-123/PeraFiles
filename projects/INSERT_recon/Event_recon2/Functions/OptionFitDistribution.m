function [options]=OptionFitDistribution(fit_method)

switch(fit_method)
    case 1 %2D Gaussian
        %opzioni funzione di fitting gaussian 2D
        options.a0=0.2; options.b0=0.2; options.c0=0.2; options.offset0=0;
        options.a_inf=0; options.b_inf=0; options.c_inf=0; options.offset_inf=-1000;
        options.a_sup=100; options.b_sup=100; options.c_sup=100; options.offset_sup=1000;

    case 2 %B-Splines
        %opzioni funzione di fitting B-spline
        options=0;
end

end


