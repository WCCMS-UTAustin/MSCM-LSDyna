      subroutine invert3x3(A, Ainv)
      implicit none
      double precision A(3,3), Ainv(3,3), det

c--- Compute determinant
      det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) -
     &      A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) +
     &      A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

c--- Prevent division by zero
      if (abs(det) .lt. 1.0d-12) then
         det = 1.0d-12
      endif

      det = 1.0d0 / det

c--- Compute inverse using cofactor method
      Ainv(1,1) =  det*(A(2,2)*A(3,3) - A(2,3)*A(3,2))
      Ainv(1,2) = -det*(A(1,2)*A(3,3) - A(1,3)*A(3,2))
      Ainv(1,3) =  det*(A(1,2)*A(2,3) - A(1,3)*A(2,2))

      Ainv(2,1) = -det*(A(2,1)*A(3,3) - A(2,3)*A(3,1))
      Ainv(2,2) =  det*(A(1,1)*A(3,3) - A(1,3)*A(3,1))
      Ainv(2,3) = -det*(A(1,1)*A(2,3) - A(1,3)*A(2,1))

      Ainv(3,1) =  det*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
      Ainv(3,2) = -det*(A(1,1)*A(3,2) - A(1,2)*A(3,1))
      Ainv(3,3) =  det*(A(1,1)*A(2,2) - A(1,2)*A(2,1))

      return
      end

