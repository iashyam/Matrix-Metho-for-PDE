program matrixMethod
    implicit NONE
    real::R,T,h,a,V0,m
    integer::l
    complex::r0,t0, DET, k1,k2,k3,e
    complex, dimension(2,2):: Mat,d21,d32,pMatrix
    h=1.
    a=1.
    v0=5.
    m=1.

    open(1, file="data.txt")
        do l=1,100
            e = cmplx(l/1000,0)
            if (E/=(0,0) .and. E/=V0) then
                k1 = csqrt(cmplx(2,0)*m*e/(h**2))
                k2 = csqrt(cmplx(2,0)/(h**2))
                k3 = k1

                !scattering matrix for k1,k2
                d21(1,1) = 0.5*(1+k1/k2)
                d21(1,2) = 0.5*(1-k1/k2)
                d21(2,1) = 0.5*(1-k1/k2)
                d21(2,2) = 0.5*(1+k1/k2)
            
                !Scattering matrix for k2, k3
                d32(1,2) = 0.5*(1-k2/k3)
                d32(2,1) = 0.5*(1-k2/k3)
                d32(2,2) = 0.5*(1+k2/k3)
                d32(1,1) = 0.5*(1+k2/k3)
            

                pMatrix(1,1) = cexp(cmplx(a,0)*cmplx(0,1)*k2)
                pMatrix(1,2) = cexp(cmplx(-a,0)*cmplx(0,1)*k2)
                pMatrix(2,1) = 0
                pMatrix(2,2) = 0

                Mat = matmul(d32, matmul(pMatrix, d21))

                r0 = -(Mat(2,1)/Mat(2,2))
                t0 = DET(Mat)/Mat(1,1)

                R = (cabs(r0))**2
                T = (cabs(t0))**2

                write(1,*) e, R, T, R+T
            endif
        end do
    close(1)
end program matrixMethod


complex function DET(M)
    implicit none
    complex,dimension(2,2),intent(in)::M
    DET = M(1,1)*M(2,2)-M(1,2)*M(2,1)
    return
end function
