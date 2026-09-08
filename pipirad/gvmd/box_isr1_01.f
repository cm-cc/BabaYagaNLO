        elmat1 = (D0aex12x45x13epem1m2p*ep3*2D0*
     &    (ame**8*ampi**2*(kp1 + kp2)*(ampi**2 - p3p4*1D0) + 
     &      p1p2*p1p3*(p1p3*2D0*
     &          (kp1**2*p2p3*(p2p4 - p2p3*1D0) + 
     &            kp2*(kp4*p1p2*p2p3 + 
     &               p1p2*
     &                (p1p4*p2p3 + p1p3*p2p4 - kp3*p2p4*1D0 - 
     &                  p1p2*p3p4*1D0 - p1p3*p2p3*2D0) + 
     &               kp2*
     &                (p1p2*p3p4 + p1p3*p2p3*2D0 - 
     &                  p1p4*p2p3*2D0)) + 
     &            kp1*p2p3*
     &             (kp3*p1p2 - kp4*p1p2*1D0 + 
     &               kp2*
     &                (p1p3 - p1p4*1D0 - p2p3*2D0 + p2p4*2D0))
     &            ) + ampi**2*p1p2*
     &          (kp1*p1p2*(kp4 - kp3*1D0) + 
     &            kp1**2*(p2p3 - p2p4*1D0) + 
     &            kp1*kp2*
     &             (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0) + 
     &            kp2*p1p2*(kp3 - kp4*1D0 + p1p3*2D0) + 
     &            kp2**2*(p1p4 - p1p3*3D0))) + 
     &      ame**6*(p1p3*p3p4*
     &          (kp2*(p1p3 + p2p3) + kp1*(p1p3 - p2p3*1D0)) + 
     &         ampi**4*((kp1 + kp2)*p1p2 + 
     &            kp1*(kp1 - kp2*3D0)) - 
     &         ampi**2*1D0*
     &          (kp1**2*p3p4 + 
     &            kp2*(p1p2*p3p4 + 
     &               (p1p3 + p2p3)*(p1p3 + p2p3 - p2p4*1D0) + 
     &               kp3*(p1p4 + p1p3*2D0) - kp4*p1p3*3D0) - 
     &            kp1*1D0*
     &             (p1p3**2 + p2p3*p2p4 - p2p3**2*1D0 - 
     &               p1p2*p3p4*1D0 - p1p3*p1p4*2D0 - 
     &               p1p3*p2p3*2D0 + 
     &               kp3*
     &                (p1p4 + p1p3*2D0 + p2p3*2D0 - p2p4*2D0)
     &                - kp4*p1p3*3D0 + p1p3*p2p4*3D0 + 
     &               kp2*p3p4*3D0))) + 
     &      ame**4*(ampi**4*p1p2*
     &          (kp2**2 - (kp1 + kp2)*p1p2*1D0 + kp1**2*3D0 - 
     &            kp1*kp2*4D0) + 
     &         p1p3*(kp1**2*p3p4*(p1p3 + p2p3*2D0) + 
     &            kp2*(p1p4*p2p3**2 - p1p3*p1p4*p2p3*1D0 - 
     &               p1p3**2*p2p4*1D0 - p1p3*p2p3*p2p4*1D0 - 
     &               kp2*p2p3*p3p4*1D0 + p1p3**2*p2p3*2D0 + 
     &               p2p3**3*2D0 - p2p3**2*p2p4*2D0 + 
     &               p1p2*p2p3*p3p4*2D0 + 
     &               kp3*p1p3*(p1p4 + p1p3*2D0) + 
     &               kp3*p2p3*(p2p4 + p2p3*2D0) - 
     &               kp4*(p1p3**2 + p2p3**2)*3D0) - 
     &            kp1*1D0*
     &             (p1p3*p2p3*p2p4 - p1p3*p1p4*p2p3*1D0 + 
     &               p1p3**3*2D0 - p1p3**2*p1p4*2D0 - 
     &               p1p3**2*p2p3*2D0 + p1p3*p2p3**2*2D0 - 
     &               p2p3**3*2D0 + p2p3**2*p2p4*2D0 - 
     &               p1p4*p2p3**2*3D0 - 
     &               kp4*(p1p3**2 + p2p3**2)*3D0 + 
     &               p1p3**2*p2p4*3D0 + 
     &               p3p4*
     &                (kp2*p2p3 - p1p2*p1p3*2D0 + 
     &                  p1p2*p2p3*2D0 + kp2*p1p3*3D0) + 
     &               kp3*
     &                (p1p3**2*2D0 + p1p3*(p1p4 - p2p4*2D0) + 
     &                  p2p3*
     &                   (-(p2p4*1D0) + p1p4*2D0 + p2p3*4D0)))
     &            ) - ampi**2*1D0*
     &          (kp2*(-(kp4*p1p2*p1p3*1D0) + 
     &               kp2*
     &                (p1p2*p3p4 - p1p4*p2p3*1D0 - 
     &                  p1p3**2*2D0 + p1p3*p1p4*2D0) + 
     &               p1p2*
     &                (-(p1p2*p3p4*1D0) + 
     &                  (kp3 + p2p3)*
     &                   (p1p4 + p2p3 - p2p4*1D0) + 
     &                  p1p3**2*3D0 + 
     &                  p1p3*
     &                   (-(p2p3*2D0) - p1p4*3D0 + p2p4*3D0)))
     &              + kp1**2*
     &             (p1p3**2 + p2p3**2 - 
     &               p2p4*1D0*(p2p3 + p1p3*2D0) + 
     &               p1p2*p3p4*3D0 + p1p3*p2p3*4D0) + 
     &            kp1*(kp4*p1p2*p1p3 + p1p2*p1p3**2 + 
     &               p1p2*p1p3*p1p4 + kp2*p1p4*p2p3 + 
     &               p1p2*p1p4*p2p3 + p1p2*p2p3**2 + 
     &               kp2*p2p3*p2p4 - kp2*p2p3**2*1D0 - 
     &               p1p2*p1p3*p2p4*1D0 - 
     &               p1p2*p2p3*p2p4*1D0 + kp2*p1p3*p1p4*2D0 - 
     &               p1p2*p1p3*p2p3*2D0 - kp2*p1p3*p2p4*2D0 + 
     &               kp3*p1p2*
     &                (p1p3*2D0 - 
     &                  (p1p4 + p2p3 - p2p4*1D0)*3D0) - 
     &               p1p2*p3p4*1D0*(p1p2 + kp2*4D0) - 
     &               kp2*p1p3**2*5D0))) + 
     &      ame**2*(-(ampi**4*p1p2**2*1D0*
     &            (kp1*(kp2 + p1p2) + kp2*(p1p2 - kp2*1D0) - 
     &              kp1**2*2D0)) + 
     &         ampi**2*p1p2*
     &          (-(kp1**2*1D0*
     &               (p1p3**2 + p2p3**2 - 
     &                 (p1p3 + p2p3)*p2p4*1D0 + 
     &                 p1p2*p3p4*2D0 + p1p3*p2p3*3D0)) - 
     &            kp2*1D0*
     &             (kp2*p1p3**2 + kp2*p1p3*p1p4 + 
     &               p1p2*p1p4*p2p3 - kp2*p1p4*p2p3*1D0 + 
     &               p1p2*p3p4*(kp2 - p1p2*1D0) + 
     &               kp4*p1p2*p1p3*3D0 - p1p2*p1p3*p1p4*3D0 + 
     &               kp3*p1p2*(p2p3 - p2p4*1D0 - p1p3*3D0) - 
     &               p1p2*p1p3*p2p3*4D0 + p1p2*p1p3*p2p4*4D0)
     &             + kp1*
     &             (p1p2*p1p3*p1p4 + kp2*p1p3*p2p3 + 
     &               kp2*p2p3**2 + kp2*p1p3*p2p4 + 
     &               p1p2*(kp2 + p1p2)*p3p4 - 
     &               kp2*p1p3*p1p4*1D0 - kp2*p1p4*p2p3*1D0 - 
     &               p1p2*p1p4*p2p3*1D0 - kp2*p2p3*p2p4*1D0 - 
     &               p1p2*p1p3**2*2D0 - p1p2*p1p3*p2p4*2D0 + 
     &               kp4*p1p2*p1p3*3D0 + kp2*p1p3**2*4D0 + 
     &               p1p2*p1p3*p2p3*4D0 + 
     &               kp3*p1p2*
     &                (p2p3 - p2p4*1D0 + p1p4*2D0 - p1p3*5D0))
     &            ) + p1p3*
     &          (kp1**2*((p1p3**2 + p1p3*p2p3 + p2p3**2)*
     &                (p2p3 - p2p4*1D0)*2D0 + 
     &               p1p2*p3p4*(p1p3 + p2p3*2D0)) + 
     &            kp2*(-(kp2*(p1p3**2 + p2p3**2)*
     &                  (p1p3 - p1p4*1D0)*2D0) + 
     &               kp2*p1p2*p3p4*(-(p2p3*1D0) + p1p3*2D0) + 
     &               p1p2*
     &                (p2p3*
     &                   (kp3*p2p4 + p1p2*p3p4 - 
     &                     p1p4*p2p3*1D0) + p1p3**3*2D0 + 
     &                  p1p3**2*(-(p1p4*2D0) + p2p4*3D0) + 
     &                  p1p3*
     &                   (-(p1p2*p3p4*3D0) + 
     &                     p2p3*
     &                     (-(p2p3*2D0) - p1p4*3D0 + p2p4*3D0)
     &                      + kp3*(p1p4 - p2p4*2D0 - p2p3*4D0)
     &                     )) - 
     &               kp4*p1p2*1D0*
     &                (p1p3**2 + p2p3**2 - p1p3*p2p3*6D0)) + 
     &            kp1*(-(kp2*1D0*
     &                  (p1p2*p2p3*p3p4 + p1p3**3*2D0 + 
     &                    p1p3*p2p3*(p2p3 - p1p4*1D0)*2D0 + 
     &                    p1p3**2*(p2p4 - p1p4*1D0)*2D0 - 
     &                    p2p3**2*(p1p4 + p2p3 - p2p4*1D0)*
     &                     2D0 + p1p2*p1p3*p3p4*3D0)) + 
     &               kp4*p1p2*
     &                (p1p3**2 + p2p3**2 - p1p3*p2p3*6D0) + 
     &               p1p2*
     &                (p1p2*p1p3*p3p4 + 
     &                  p2p3*(p1p4*p2p3 - p1p2*p3p4*1D0) + 
     &                  kp3*p2p3*(p2p4 - (p1p4 + p2p3)*2D0) + 
     &                  p1p3*p2p3*
     &                   (-(p1p4*3D0) + p2p4*3D0 - p2p3*4D0)
     &                   + p1p3**2*(-(p2p4*1D0) + p2p3*4D0) + 
     &                  kp3*p1p3*
     &                   (-(p1p4*1D0) + p2p4*2D0 + p2p3*6D0)))
     &            ))))/
     &  (kp1*kp2*(ampi**2 + p3p4)*
     &    (ame**4*ampi**2 - ame**2*(p1p3**2 + p2p3**2)*1D0 + 
     &      p1p2*(-(ampi**2*p1p2*1D0) + p1p3*p2p3*2D0)))
                  
