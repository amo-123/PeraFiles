function [Frame] = mask_frame( Frame )
%Neighbourhood mask, 21 closest neighbours

    nx = 6;  ny = 12;
    nv = size(Frame,1);
    [~,ii_max] = max(Frame,[],2);
    ix = mod( ii_max-1, nx ) + 1;
    iy = fix( (ii_max-1) / nx ) + 1;
    msk = zeros([nx,ny,nv]);
    ix1=ix-2; ix1(ix1<1)=1; ix2=ix+2; ix2(ix2>nx)=nx;
    iy1=iy-1; iy1(iy1<1)=1; iy2=iy+1; iy2(iy2>ny)=ny;
    jx1=ix-1; jx1(jx1<1)=1; jx2=ix+1; jx2(jx2>nx)=nx;
    jy1=iy-2; jy1(jy1<1)=1; jy2=iy+2; jy2(jy2>ny)=ny;
    for iv=1:nv
        msk( ix1(iv):ix2(iv), iy1(iv):iy2(iv), iv )=1;
        msk( jx1(iv):jx2(iv), jy1(iv):jy2(iv), iv )=1;
    end
    msk = reshape( msk, [nx*ny,nv] )';
    Frame = Frame .* msk;
    
end
