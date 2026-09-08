        elmat1 = -((D0aex12x45x13epem1m2p*ep3*1D0*
     &      (kp1**2*(ame**2*
     &            (ampi**2*p2p3**2 - ampi**2*p2p3*p2p4*1D0 + 
     &              ame**2*ampi**2*(p3p4 - ampi**2*1D0)) + 
     &           p1p3**2*
     &            (ame**2*(ampi**2 - p3p4*1D0) - 
     &              p2p3**2*2D0 + p2p3*p2p4*2D0) + 
     &           ame**2*p1p3*p2p3*(ampi**2*2D0 - p3p4*2D0) + 
     &           p1p2*(p1p3*
     &               (ampi**2*p2p3 - ampi**2*p2p4*1D0) - 
     &              ame**2*ampi**2*(ampi**2 - p3p4*1D0)*2D0))
     &         + kp1*(ame**2*ampi**4*p1p2**2 + 
     &           ampi**2*kp4*p1p2**2*p1p3 + 
     &           ame**4*ampi**2*p1p3**2 + 
     &           ame**2*kp4*p1p3**3 + 
     &           ame**2*ampi**2*p1p2*p1p3*p1p4 + 
     &           ame**2*ampi**2*p1p2*p1p4*p2p3 + 
     &           ame**4*ampi**2*p2p3**2 + 
     &           ame**2*kp4*p1p3*p2p3**2 + 
     &           ame**2*p1p3**3*p2p4 + 
     &           ame**2*p1p3**2*p2p3*p2p4 + 
     &           ame**6*ampi**2*p3p4 + 
     &           ame**4*p1p3*p2p3*p3p4 + 
     &           ame**2*p1p2*p1p3*p2p3*p3p4 - 
     &           ame**6*ampi**4*1D0 - 
     &           ame**4*ampi**2*kp4*p1p3*1D0 - 
     &           ame**2*p1p3**2*p1p4*p2p3*1D0 - 
     &           ame**2*p1p3*p1p4*p2p3**2*1D0 - 
     &           ame**4*ampi**2*p1p3*p2p4*1D0 - 
     &           ame**4*ampi**2*p2p3*p2p4*1D0 - 
     &           ame**2*ampi**2*p1p2**2*p3p4*1D0 - 
     &           ame**4*p1p3**2*p3p4*1D0 - 
     &           ame**2*p1p2*p1p3**2*p3p4*1D0 - 
     &           ame**2*ampi**2*p1p2*p1p3*p2p3*2D0 - 
     &           kp4*p1p2*p1p3**2*p2p3*2D0 - 
     &           kp3*1D0*
     &            (ampi**2*p1p2**2*p1p3 + 
     &              p1p2*
     &               (-(ame**2*ampi**2*p1p3*2D0) - 
     &                 p1p3**2*p2p3*2D0 + 
     &                 ame**2*
     &                  (ampi**2*p2p3 - ampi**2*p2p4*1D0 + 
     &                    ampi**2*p1p4*2D0)) + 
     &              ame**2*
     &               (ame**2*ampi**2*p1p4 + p1p3**3*2D0 + 
     &                 ame**2*ampi**2*p2p3*2D0 - 
     &                 ame**2*ampi**2*p2p4*2D0 - 
     &                 p1p3*1D0*
     &                  (-(p2p3*p2p4*1D0) + 
     &                    ame**2*ampi**2*2D0 + p1p4*p2p3*2D0)
     &                  + p1p3**2*(-(p1p4*1D0) + p2p4*2D0)))
     &            + kp2*(p1p2*
     &               (-(ampi**2*p1p3**2*1D0) + 
     &                 p1p3*
     &                  (ampi**2*p1p4 + ampi**2*p2p3 - 
     &                    ampi**2*p2p4*1D0) + 
     &                 ame**2*ampi**2*(ampi**2 - p3p4*1D0)) + 
     &              p1p3**3*p2p3*2D0 + 
     &              ame**2*p1p3*p2p3*(p3p4 - ampi**2*2D0) + 
     &              ame**2*
     &               (ampi**2*p1p4*p2p3 + ampi**2*p2p3*p2p4 - 
     &                 ampi**2*p2p3**2*1D0 + 
     &                 ame**2*ampi**2*(ampi**2 - p3p4*1D0)*3D0
     &                 ) + 
     &              p1p3**2*
     &               (-(p1p4*p2p3*2D0) + 
     &                 ame**2*(-(ampi**2*3D0) + p3p4*3D0))))
     &         + kp2*(ame**2*ampi**4*p1p2**2 + 
     &           ame**4*ampi**2*kp4*p1p3 + 
     &           ame**2*ampi**2*p1p2*p1p4*p2p3 + 
     &           ame**2*p1p3**2*p1p4*p2p3 + 
     &           ame**4*ampi**2*p2p3**2 + 
     &           ame**2*p1p3*p1p4*p2p3**2 + 
     &           ame**4*ampi**2*p1p3*p2p4 + 
     &           ame**2*p1p3**2*p2p3*p2p4 + 
     &           ame**6*ampi**2*p3p4 + 
     &           ame**2*p1p2*p1p3**2*p3p4 - 
     &           ame**6*ampi**4*1D0 - 
     &           ampi**2*kp4*p1p2**2*p1p3*1D0 - 
     &           ame**2*kp4*p1p3**3*1D0 - 
     &           ame**2*ampi**2*p1p2*p1p3*p1p4*1D0 - 
     &           ame**2*kp4*p1p3*p2p3**2*1D0 - 
     &           ame**2*p1p3**3*p2p4*1D0 - 
     &           ame**4*ampi**2*p2p3*p2p4*1D0 - 
     &           ame**2*ampi**2*p1p2**2*p3p4*1D0 - 
     &           ame**4*p1p3**2*p3p4*1D0 - 
     &           ame**4*p1p3*p2p3*p3p4*1D0 - 
     &           ame**2*p1p2*p1p3*p2p3*p3p4*1D0 - 
     &           ampi**2*p1p2**2*p1p3**2*2D0 - 
     &           ame**2*p1p3**4*2D0 - 
     &           ame**4*ampi**2*p1p3*p1p4*2D0 + 
     &           ame**2*p1p3**3*p1p4*2D0 - 
     &           ame**2*ampi**2*p1p2*p1p3*p2p3*2D0 + 
     &           kp4*p1p2*p1p3**2*p2p3*2D0 - 
     &           p1p2*p1p3**2*p1p4*p2p3*2D0 - 
     &           ame**2*p1p3**2*p2p3**2*2D0 + 
     &           ame**2*ampi**2*p1p2*p1p3*p2p4*2D0 - 
     &           p1p2*p1p3**3*p2p4*2D0 + 
     &           p1p2**2*p1p3**2*p3p4*2D0 + 
     &           kp3*(ampi**2*p1p2**2*p1p3 + 
     &              p1p2*
     &               (ame**2*
     &                  (ampi**2*p2p3 - ampi**2*p2p4*1D0) + 
     &                 p1p3**2*2D0*(p2p4 - p2p3*2D0)) + 
     &              ame**2*
     &               (ame**2*ampi**2*p1p4 - 
     &                 p1p3**2*p1p4*1D0 + p1p3**3*2D0 - 
     &                 p1p3*1D0*
     &                  (p2p3*p2p4 + ame**2*ampi**2*2D0 - 
     &                    p2p3**2*2D0))) + 
     &           kp2*(ame**2*
     &               (p1p3*p2p3*p3p4 - ampi**2*p1p4*p2p3*1D0)
     &               + p1p2*
     &               (ampi**2*p1p3*p1p4 - 
     &                 ame**2*ampi**2*1D0*
     &                  (ampi**2 - p3p4*1D0) + 
     &                 p1p3**2*(ampi**2 - p3p4*2D0))) + 
     &           ame**4*ampi**2*p1p3**2*3D0 + 
     &           p1p2*p1p3**3*p2p3*4D0)))/
     &    (kp1*kp2*(ampi**2 + p3p4)*
     &      (ampi**2*p1p2**2 + 
     &        ame**2*(p1p3**2 + p2p3**2 - 
     &           ame**2*ampi**2*1D0) - p1p2*p1p3*p2p3*2D0)))
                  
