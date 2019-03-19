function outputs = FrangiFilter3D_wrapper(I,options)
    [Iout,whatScale,Voutx,Vouty,Voutz]=FrangiFilter3D(I,options);
    outputs{1} = Iout;
    outputs{2} = whatScale;
end