function dt = CFLcond(aAdap, alphaAdap, drDense, fg)
dt = fg*drDense./(max(aAdap./alphaAdap));
end