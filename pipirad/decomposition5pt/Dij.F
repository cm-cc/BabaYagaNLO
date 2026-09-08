      bc0 = C0BCD + D0m*m1q
               
      bc1 = (C0BCD*k1k1 + C1ACD*k1k2 - C1BCD*(-k1k1 + k1k2) + C2ACD*k1k3
     >     - C2BCD*(-k1k1 + k1k3) + (D1m*k1k1 + D2m*k1k2 + D3m*k1k3)*
     >     (-k1k1 - m1q + m2q))/2.q0
      
      bc2 = (C0BCD*k1k1 + C1ABD*k1k1 - C1BCD*(-k1k1 + k1k2) + C2ABD*k1k3
     >     - C2BCD*(-k1k1 + k1k3) + (D1m*k1k1 + D2m*k1k2 + D3m*k1k3)*
     >     (-k2k2 - m1q + m3q))/2.q0

!      print*,'B2',bc2
      bc2 = (C0BCD*k1k2 + C1ACD*k2k2 - C1BCD*(-k1k2 + k2k2) + C2ACD*k2k3
     >     - C2BCD*(-k1k2 + k2k3) + (D1m*k1k2 + D2m*k2k2 + D3m*k2k3)*
     >     (-k1k1 - m1q + m2q))/2.q0 +bc2

      bc2 = bc2/2.q0
      
      
      bc3 = (C0BCD*k1k1 + C1ABC*k1k1 + C2ABC*k1k2 - C1BCD*(-k1k1 + k1k2)
     >     - C2BCD*(-k1k1 + k1k3) + (D1m*k1k1 + D2m*k1k2 + D3m*k1k3)
     >     *(-k3k3 - m1q + m4q))/2.q0

      bc3 = (C0BCD*k1k3 + C1ACD*k2k3 - C1BCD*(-k1k3 + k2k3) + C2ACD*k3k3
     >     - C2BCD*(-k1k3 + k3k3) + (D1m*k1k3 + D2m*k2k3 + D3m*k3k3)*
     >     (-k1k1 - m1q + m2q))/2.q0 + bc3

      bc3 = bc3/2.q0
      
      bc4 = (C0BCD*k1k2 + C1ABD*k1k2 - C1BCD*(-k1k2 + k2k2) + C2ABD*k2k3
     >     - C2BCD*(-k1k2 + k2k3) + (D1m*k1k2 + D2m*k2k2 + D3m*k2k3)*
     >     (-k2k2 - m1q + m3q))/2.q0

      
      bc5 = (C0BCD*k1k3 + C1ABD*k1k3 - C1BCD*(-k1k3 + k2k3) + C2ABD*k3k3
     >     - C2BCD*(-k1k3 + k3k3) + (D1m*k1k3 + D2m*k2k3 + D3m*k3k3)*
     >     (-k2k2 - m1q + m3q))/2.q0

      
      bc5 = (C0BCD*k1k2 + C1ABC*k1k2 + C2ABC*k2k2 - C1BCD*(-k1k2 + k2k2)
     >     - C2BCD*(-k1k2 + k2k3) + (D1m*k1k2 + D2m*k2k2 + D3m*k2k3)*
     >     (-k3k3 - m1q + m4q))/2.q0 +bc5

      bc5 =bc5/2.q0
      
      bc6 = (C0BCD*k1k3 + C1ABC*k1k3 + C2ABC*k2k3 - C1BCD*(-k1k3 + k2k3)
     >     - C2BCD*(-k1k3 + k3k3) + (D1m*k1k3 + D2m*k2k3 + D3m*k3k3)*
     >     (-k3k3 - m1q + m4q))/2.q0

      
      Dijm(1) = (bc6*k1k1*k2k2 + bc0*k1k3**2*k2k2 + 
     &    bc0*k1k1*k2k3**2 - bc6*k1k2**2*1.q0 - 
     &    bc4*k1k3**2*1.q0 - bc1*k2k3**2*1.q0 + 
     &    bc5*k1k2*k1k3*2.q0 - bc3*k1k3*k2k2*2.q0 - 
     &    bc5*k1k1*k2k3*2.q0 + bc3*k1k2*k2k3*2.q0 + 
     &    bc2*k1k3*k2k3*2.q0 - bc0*k1k2*k1k3*k2k3*2.q0 + 
     &    k3k3*(bc4*k1k1 + bc0*k1k2**2 + bc1*k2k2 - 
     &       bc0*k1k1*k2k2*1.q0 - bc2*k1k2*2.q0))/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &    k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0)
        Dijm(2) = (bc6*(k1k3**2*k2k2**2 + 
     &       k2k2*k3k3*(k1k1*k2k2 - k1k2**2*1.q0) - 
     &       k1k2*k1k3*k2k2*k2k3*2.q0 + 
     &       k2k3**2*(-(k1k1*k2k2*1.q0) + k1k2**2*2.q0)) + 
     &    k3k3**2*(bc4*(k1k2**2 + k1k1*k2k2) + 
     &       k2k2*(bc0*k1k2**2 - bc0*k1k1*k2k2*1.q0 + 
     &          bc1*k2k2*2.q0 - bc2*k1k2*4.q0)) + 
     &    k2k3*(-(bc5*k1k3**2*k2k2*2.q0) + 
     &       k2k3**3*(-(bc0*k1k1*1.q0) + bc1*2.q0) + 
     &       k2k3**2*2.q0*
     &        (bc5*k1k1 + bc0*k1k2*k1k3 - bc3*k1k2*2.q0 - 
     &          bc2*k1k3*2.q0) + 
     &       k1k3*k2k3*(-(bc0*k1k3*k2k2*1.q0) + bc4*k1k3*2.q0 + 
     &          bc3*k2k2*4.q0)) - 
     &    k3k3*1.q0*(k2k3*2.q0*
     &        (k1k2*(bc5*k1k2 + bc4*k1k3) + 
     &          k2k2*(bc5*k1k1 + bc0*k1k2*k1k3 - 
     &             bc3*k1k2*2.q0 - bc2*k1k3*2.q0)) + 
     &       k2k3**2*(bc4*k1k1 + bc0*k1k2**2 - 
     &          bc0*k1k1*k2k2*2.q0 - bc2*k1k2*4.q0 + 
     &          bc1*k2k2*4.q0) + 
     &       k1k3*k2k2*(bc4*k1k3 - bc0*k1k3*k2k2*1.q0 - 
     &          bc5*k1k2*4.q0 + bc3*k2k2*4.q0)))/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &     k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0
     &     )**2
        Dijm(5) = (k1k3**3*k2k3*2.q0*(bc0*k1k2 - bc2*2.q0) - 
     &    k1k3*k2k3**2*1.q0*
     &     (bc0*k1k1*k1k3 + bc3*k1k1*2.q0 - bc1*k1k3*2.q0) + 
     &    k1k3**3*(-(bc0*k1k3*k2k2*1.q0) + bc4*k1k3*2.q0 + 
     &       bc3*k2k2*2.q0) + 
     &    bc6*(k1k1*(k1k1*(k2k3**2 + k2k2*k3k3) - 
     &          k1k3**2*k2k2*1.q0) - k1k1*k1k2*k1k3*k2k3*2.q0 + 
     &       k1k2**2*(-(k1k1*k3k3*1.q0) + k1k3**2*2.q0)) - 
     &    bc5*(k1k2*k1k3 - k1k1*k2k3*1.q0)*
     &     (k1k3**2 - k1k1*k3k3*1.q0)*4.q0 + 
     &    k3k3**2*((bc1 + bc0*k1k1)*k1k2**2 + 
     &       k1k1*k2k2*(bc1 - bc0*k1k1*1.q0) + 
     &       bc4*k1k1**2*2.q0 - bc2*k1k1*k1k2*4.q0) - 
     &    k3k3*1.q0*(k1k1*k2k3**2*(bc1 - bc0*k1k1*1.q0) + 
     &       bc3*k1k3*(k1k2**2 + k1k1*k2k2)*2.q0 + 
     &       k1k3*k2k3*2.q0*
     &        ((bc1 + bc0*k1k1)*k1k2 - bc2*k1k1*2.q0) - 
     &       bc3*k1k1*k1k2*k2k3*4.q0 + 
     &       k1k3**2*(bc0*k1k2**2 + bc1*k2k2 - 
     &          bc0*k1k1*k2k2*2.q0 + bc4*k1k1*4.q0 - 
     &          bc2*k1k2*4.q0)))/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &     k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0
     &     )**2
        Dijm(7) = ((bc4*k1k1 + bc0*k1k2**2 + bc1*k2k2 - 
     &       bc0*k1k1*k2k2*1.q0 - bc2*k1k2*2.q0)*
     &     (k1k1*k2k2*k3k3 - k1k3**2*k2k2*1.q0 - 
     &       k1k1*k2k3**2*1.q0 - k1k2**2*k3k3*1.q0 + 
     &       k1k2*k1k3*k2k3*2.q0) + 
     &    2.q0*(bc6*(k1k2**2 - k1k1*k2k2*1.q0)**2 + 
     &       k2k3*(k1k2*(bc5*k1k1*k1k2 + bc2*k1k2*k1k3 - 
     &             bc3*k1k2**2*1.q0 - bc4*k1k1*k1k3*1.q0) + 
     &          k2k2*(bc3*k1k1*k1k2 + bc2*k1k1*k1k3 - 
     &             bc5*k1k1**2*1.q0 - bc1*k1k2*k1k3*1.q0))*2.q0
     &        + k2k3**2*(bc4*k1k1**2 + 
     &          k1k2*(bc1*k1k2 - bc2*k1k1*2.q0)) + 
     &       k1k3*(k1k2*k2k2*
     &           (bc5*k1k1 + bc3*k1k2 - bc2*k1k3*1.q0)*2.q0 + 
     &          k2k2**2*(bc1*k1k3 - bc3*k1k1*2.q0) + 
     &          k1k2**2*(bc4*k1k3 - bc5*k1k2*2.q0))))/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &     k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0
     &     )**2
        Dijm(3) = (bc6*(k1k2**3*k3k3 - 
     &       k1k2*(k1k3**2*k2k2 + k1k1*k2k3**2 + 
     &          k1k1*k2k2*k3k3)*1.q0 + k1k1*k1k3*k2k2*k2k3*2.q0)
     &      + k3k3**2*(bc2*k1k1*k2k2 + bc0*k1k1*k1k2*k2k2 - 
     &       bc0*k1k2**3*1.q0 - bc4*k1k1*k1k2*2.q0 - 
     &       bc1*k1k2*k2k2*2.q0 + bc2*k1k2**2*3.q0) + 
     &    k2k3*(k2k3**2*(bc3*k1k1 + bc0*k1k1*k1k3 - 
     &          bc1*k1k3*2.q0) + 
     &       k1k3*k2k3*2.q0*
     &        (bc3*k1k2 - bc0*k1k2*k1k3*1.q0 + bc2*k1k3*2.q0) + 
     &       k1k3**2*(bc0*k1k3*k2k2 - bc4*k1k3*2.q0 - 
     &          bc3*k2k2*3.q0)) + 
     &    bc5*(k1k3**3*k2k2 + k1k2*k1k3**2*k2k3*2.q0 - 
     &       k1k3*1.q0*(k1k1*k2k2*k3k3 + k1k1*k2k3**2*3.q0 + 
     &          k1k2**2*k3k3*3.q0) + k1k1*k1k2*k2k3*k3k3*4.q0)
     &     + k3k3*(bc4*k1k3*(k1k2*k1k3 + k1k1*k2k3)*2.q0 - 
     &       k2k3**2*1.q0*
     &        (bc2*k1k1 + bc0*k1k1*k1k2 - bc1*k1k2*2.q0) - 
     &       k2k3*1.q0*(k2k2*
     &           (bc3*k1k1 + bc0*k1k1*k1k3 - bc1*k1k3*2.q0) + 
     &          k1k2*(bc3*k1k2 - bc0*k1k2*k1k3*1.q0 + 
     &             bc2*k1k3*2.q0)*3.q0) - 
     &       k1k3*k2k2*1.q0*
     &        (bc2*k1k3 + bc0*k1k2*k1k3 - bc3*k1k2*4.q0)))/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &     k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0
     &     )**2
        Dijm(4) = (bc4*k1k3**3*k2k2 + bc2*k1k1*k2k3**3 + 
     &    bc0*k1k1*k1k2*k2k3**3 - bc0*k1k3**3*k2k2**2*1.q0 - 
     &    bc4*k1k1*k1k3*k2k3**2*1.q0 - 
     &    bc3*k1k1*k2k2*k2k3**2*1.q0 - 
     &    bc0*k1k1*k1k3*k2k2*k2k3**2*1.q0 + 
     &    bc5*k1k2**2*k1k3*k2k3*2.q0 + 
     &    bc2*k1k2*k1k3*k2k3**2*2.q0 - 
     &    bc0*k1k2**2*k1k3*k2k3**2*2.q0 + 
     &    bc1*k1k3*k2k2*k2k3**2*2.q0 - bc1*k1k2*k2k3**3*2.q0 - 
     &    bc6*(k1k2**2 - k1k1*k2k2*1.q0)*
     &     (k1k2*k2k3 - k1k3*k2k2*1.q0)*2.q0 - 
     &    bc5*k1k2*k1k3**2*k2k2*3.q0 + 
     &    bc3*k1k3**2*k2k2**2*3.q0 - 
     &    bc2*k1k3**2*k2k2*k2k3*3.q0 + 
     &    bc0*k1k2*k1k3**2*k2k2*k2k3*3.q0 - 
     &    bc5*k1k1*k1k2*k2k3**2*3.q0 + 
     &    bc5*k1k1*k1k3*k2k2*k2k3*4.q0 + 
     &    bc3*k1k2**2*k2k3**2*4.q0 + 
     &    k3k3*(-(bc4*k1k3*(k1k2**2 + k1k1*k2k2)*1.q0) + 
     &       bc5*(k1k2**3 - k1k1*k1k2*k2k2*1.q0) + 
     &       bc4*k1k1*k1k2*k2k3*2.q0 - 
     &       k2k3*1.q0*(bc2*k1k1*k2k2 + bc0*k1k1*k1k2*k2k2 - 
     &          bc0*k1k2**3*1.q0 - bc1*k1k2*k2k2*2.q0 + 
     &          bc2*k1k2**2*3.q0) + 
     &       k2k2*(k2k2*(bc3*k1k1 + bc0*k1k1*k1k3 - 
     &             bc1*k1k3*2.q0) - 
     &          k1k2*1.q0*
     &           (bc3*k1k2 + bc0*k1k2*k1k3 - bc2*k1k3*4.q0)))
     &     - bc3*k1k2*k1k3*k2k2*k2k3*6.q0)/
     &  (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &     k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - k1k2*k1k3*k2k3*2.q0
     &     )**2
        Dijm(6) = -((1.q0*(bc1*k1k3**2*k2k2*k2k3 + 
     &        bc0*k1k1*k1k3**2*k2k2*k2k3 + 
     &        bc0*k1k1**2*k2k3**3 - bc2*k1k3**3*k2k2*1.q0 - 
     &        bc0*k1k2*k1k3**3*k2k2*1.q0 - 
     &        bc1*k1k1*k2k3**3*1.q0 + bc4*k1k2*k1k3**3*2.q0 - 
     &        bc3*k1k2**2*k1k3*k2k3*2.q0 - 
     &        bc4*k1k1*k1k3**2*k2k3*2.q0 - 
     &        bc2*k1k2*k1k3**2*k2k3*2.q0 + 
     &        bc0*k1k2**2*k1k3**2*k2k3*2.q0 + 
     &        bc6*(k1k2**2 - k1k1*k2k2*1.q0)*
     &         (k1k2*k1k3 - k1k1*k2k3*1.q0)*2.q0 + 
     &        bc3*k1k2*k1k3**2*k2k2*3.q0 + 
     &        bc3*k1k1*k1k2*k2k3**2*3.q0 + 
     &        bc2*k1k1*k1k3*k2k3**2*3.q0 - 
     &        bc0*k1k1*k1k2*k1k3*k2k3**2*3.q0 - 
     &        bc3*k1k1*k1k3*k2k2*k2k3*4.q0 + 
     &        k3k3*(bc3*(k1k1*k1k2*k2k2 - k1k2**3*1.q0) + 
     &           k1k3*(bc2*k1k1*k2k2 + bc0*k1k1*k1k2*k2k2 - 
     &              bc0*k1k2**3*1.q0 - bc4*k1k1*k1k2*2.q0 - 
     &              bc1*k1k2*k2k2*2.q0 + bc2*k1k2**2*3.q0) + 
     &           k2k3*((bc1 + bc0*k1k1)*k1k2**2 + 
     &              k1k1*k2k2*(bc1 - bc0*k1k1*1.q0) + 
     &              bc4*k1k1**2*2.q0 - bc2*k1k1*k1k2*4.q0)) + 
     &        bc5*(k1k1*(k1k3**2*k2k2 - 
     &              k1k1*1.q0*(k2k2*k3k3 + k2k3**2*3.q0)) + 
     &           k1k2**2*(k1k1*k3k3 - k1k3**2*4.q0) + 
     &           k1k1*k1k2*k1k3*k2k3*6.q0)))/
     &    (k1k3**2*k2k2 + k1k2**2*k3k3 + 
     &       k1k1*(k2k3**2 - k2k2*k3k3*1.q0) - 
     &       k1k2*k1k3*k2k3*2.q0)**2)
