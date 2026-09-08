        elmat2 = -((A0ap*p3p4*1.25D-1*
     &       (-(x12*1D0) + (ame**2 + p1p2)*2D0)*
     &       (kp2*(kp4*p1p3 + kp3*p1p4) + 
     &         kp1*(kp4*p2p3 - kp3*p2p4*2D0)))/
     &     (ampi**2*kp1*kp2*kp3*kp4*(ame**2 + p1p2)*
     &       (ampi**2 + p3p4))) + 
     &  (A0am1*p3p4*1.25D-1*
     &     (-(x12*1D0) + (ame**2 + p1p2)*2D0)*
     &     (kp2*(kp4*p1p3 - kp3*p1p4*1D0) + 
     &       kp1*(kp4*p2p3 - kp3*p2p4*2D0)))/
     &   (ampi**2*kp1*kp2*kp3*kp4*(ame**2 + p1p2)*
     &     (ampi**2 + p3p4)) - 
     &  (C0appx34pm1p*2.5D-1*(ampi**2 - m1**2*1D0 - p3p4*1D0)*
     &     (kp1**2*(kp4*p2p3*
     &           (p2p4 + p3p4 - ampi**2*1D0 - p2p3*1D0) + 
     &          kp2*(kp3 - kp4*1D0)*(ampi**2 - p3p4*1D0) + 
     &          kp3*p2p4*
     &           (ampi**2 + p2p4 - p2p3*1D0 - p3p4*1D0)) + 
     &       kp2*(-(kp3**2*1D0*
     &             (p1p4*(ame**2 + p1p2 + p2p4 - p2p3*1D0) + 
     &               ame**2*kp4*2D0)) + 
     &          kp4*p1p3*
     &           (kp2*(ampi**2 + p1p3 - p1p4*1D0 - 
     &                p3p4*1D0) + 
     &             (p2p3 - p2p4*1D0)*
     &              (kp4 + p1p3*2D0 - p1p4*2D0) + 
     &             ame**2*(kp4 - ampi**2*2D0 + p3p4*2D0) + 
     &             p1p2*(kp4 - ampi**2*2D0 + p3p4*2D0)) + 
     &          kp3*(ame**2*kp4**2*2D0 + 
     &             kp4*(p1p4*
     &                 (p1p2 + p2p4 - ame**2*1D0 - p2p3*1D0)
     &                 + p1p3*
     &                 (ame**2 + p2p4 - p1p2*1D0 - p2p3*1D0)
     &                 + ame**2*(p2p4 - p2p3*1D0)*2D0) + 
     &             p1p4*(-(p3p4*x12*1D0) + 
     &                kp2*
     &                 (p1p3 + p3p4 - ampi**2*1D0 - p1p4*1D0)
     &                 + ame**2*ampi**2*2D0 + 
     &                ampi**2*p1p2*2D0 - 
     &                (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0)
     &             )) + kp1*
     &        (-(kp2**2*1D0*(kp3 - kp4*1D0)*
     &             (ampi**2 - p3p4*1D0)) + 
     &          kp3**2*(p2p4*
     &              (ame**2 + p1p2 + p1p4 - p1p3*1D0) + 
     &             ame**2*kp4*2D0) - 
     &          kp4*p2p3*1D0*
     &           (p3p4*x12 + 
     &             kp4*(ame**2 + p1p2 + p1p3 - p1p4*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0) + 
     &          kp3*(-(ame**2*kp4**2*2D0) + 
     &             p2p4*(p3p4*x12 - ame**2*ampi**2*2D0 - 
     &                ampi**2*p1p2*2D0 + 
     &                (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0)
     &               + kp4*
     &              (-(1D0*(ame**2 - p1p2*1D0)*
     &                   (p2p3 - p2p4*1D0)) + 
     &                p1p3*(p2p3 + p2p4 + ame**2*2D0) - 
     &                p1p4*1D0*(p2p3 + p2p4 + ame**2*2D0))) + 
     &          kp2*(kp3**2*(p1p4 - p2p4*1D0) + 
     &             kp4*(kp4*(p2p3 - p1p3*1D0) - 
     &                p2p4*1D0*(ampi**2 + p1p3 - p3p4*1D0) + 
     &                p1p4*(ampi**2 + p2p3 - p3p4*1D0)) + 
     &             kp3*(p2p3*(ampi**2 + p1p4 - p3p4*1D0) - 
     &                kp4*(p1p4 + p2p3 - p2p4*1D0)*3D0 - 
     &                p1p3*1D0*
     &                 (ampi**2 + p2p4 - p3p4*1D0 - kp4*3D0)))
     &          ))*(m1**2 + p3p4*4D0))/
     &   (kp1*kp2*kp3*kp4*(ame**2 + p1p2)*
     &     (ampi**4 - p3p4**2*1D0)) + 
     &  (B0ax34pp*2.5D-1*
     &     (kp1**2*(kp4*p2p3*
     &           (p2p4 + p3p4 - ampi**2*1D0 - p2p3*1D0) + 
     &          kp2*(kp3 - kp4*1D0)*(ampi**2 - p3p4*1D0) + 
     &          kp3*p2p4*
     &           (ampi**2 + p2p4 - p2p3*1D0 - p3p4*1D0)) + 
     &       kp2*(-(kp3**2*1D0*
     &             (p1p4*(ame**2 + p1p2 + p2p4 - p2p3*1D0) + 
     &               ame**2*kp4*2D0)) + 
     &          kp4*p1p3*
     &           (kp2*(ampi**2 + p1p3 - p1p4*1D0 - 
     &                p3p4*1D0) + 
     &             (p2p3 - p2p4*1D0)*
     &              (kp4 + p1p3*2D0 - p1p4*2D0) + 
     &             ame**2*(kp4 - ampi**2*2D0 + p3p4*2D0) + 
     &             p1p2*(kp4 - ampi**2*2D0 + p3p4*2D0)) + 
     &          kp3*(ame**2*kp4**2*2D0 + 
     &             kp4*(p1p4*
     &                 (p1p2 + p2p4 - ame**2*1D0 - p2p3*1D0)
     &                 + p1p3*
     &                 (ame**2 + p2p4 - p1p2*1D0 - p2p3*1D0)
     &                 + ame**2*(p2p4 - p2p3*1D0)*2D0) + 
     &             p1p4*(kp2*
     &                 (p1p3 + p3p4 - ampi**2*1D0 - p1p4*1D0)
     &                 - (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*
     &                 2D0 + 
     &                ame**2*(ampi**2 - p3p4*1D0)*2D0 + 
     &                p1p2*(ampi**2 - p3p4*1D0)*2D0))) + 
     &       kp1*(-(kp2**2*1D0*(kp3 - kp4*1D0)*
     &             (ampi**2 - p3p4*1D0)) + 
     &          kp3**2*(p2p4*
     &              (ame**2 + p1p2 + p1p4 - p1p3*1D0) + 
     &             ame**2*kp4*2D0) - 
     &          kp4*p2p3*1D0*
     &           (p3p4*x12 + 
     &             kp4*(ame**2 + p1p2 + p1p3 - p1p4*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0) + 
     &          kp3*(-(ame**2*kp4**2*2D0) + 
     &             p2p4*(p3p4*x12 - ame**2*ampi**2*2D0 - 
     &                ampi**2*p1p2*2D0 + 
     &                (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0)
     &               + kp4*
     &              (-(1D0*(ame**2 - p1p2*1D0)*
     &                   (p2p3 - p2p4*1D0)) + 
     &                p1p3*(p2p3 + p2p4 + ame**2*2D0) - 
     &                p1p4*1D0*(p2p3 + p2p4 + ame**2*2D0))) + 
     &          kp2*(kp3**2*(p1p4 - p2p4*1D0) + 
     &             kp4*(kp4*(p2p3 - p1p3*1D0) - 
     &                p2p4*1D0*(ampi**2 + p1p3 - p3p4*1D0) + 
     &                p1p4*(ampi**2 + p2p3 - p3p4*1D0)) + 
     &             kp3*(p2p3*(ampi**2 + p1p4 - p3p4*1D0) - 
     &                kp4*(p1p4 + p2p3 - p2p4*1D0)*3D0 - 
     &                p1p3*1D0*
     &                 (ampi**2 + p2p4 - p3p4*1D0 - kp4*3D0)))
     &          ))*(m1**2 + p3p4*4D0))/
     &   (kp1*kp2*kp3*kp4*(ame**2 + p1p2)*
     &     (ampi**4 - p3p4**2*1D0)) - 
     &  (B0appm1*1.25D-1*
     &     (ampi**2*kp1**2*
     &        (kp4*p2p3*(p2p4 + p3p4 - ampi**2*1D0 - 
     &             p2p3*1D0) + 
     &          kp2*(kp3 - kp4*1D0)*(ampi**2 - p3p4*1D0) + 
     &          kp3*p2p4*
     &           (ampi**2 + p2p4 - p2p3*1D0 - p3p4*1D0))*2D0*
     &        (m1**2 + ampi**2*2D0 + p3p4*2D0) + 
     &       kp1*(-(ampi**2*kp2**2*(kp3 - kp4*1D0)*
     &             (ampi**2 - p3p4*1D0)*2D0*
     &             (m1**2 + ampi**2*2D0 + p3p4*2D0)) + 
     &          ampi**2*kp3**2*2D0*
     &           (p2p4*(ame**2 + p1p2 + p1p4 - p1p3*1D0) + 
     &             ame**2*kp4*2D0)*
     &           (m1**2 + ampi**2*2D0 + p3p4*2D0) - 
     &          ampi**2*kp2*2D0*
     &           (m1**2 + ampi**2*2D0 + p3p4*2D0)*
     &           (kp3**2*(p2p4 - p1p4*1D0) + 
     &             kp4*(kp4*(p1p3 - p2p3*1D0) + 
     &                p2p4*(ampi**2 + p1p3 - p3p4*1D0) - 
     &                p1p4*1D0*(ampi**2 + p2p3 - p3p4*1D0)) + 
     &             kp3*(-(p2p3*1D0*
     &                   (ampi**2 + p1p4 - p3p4*1D0)) + 
     &                kp4*(p1p4 + p2p3 - p2p4*1D0)*3D0 + 
     &                p1p3*
     &                 (ampi**2 + p2p4 - p3p4*1D0 - kp4*3D0)))
     &            - kp3*2D0*
     &           (ame**2*ampi**2*kp4**2*2D0*
     &              (m1**2 + ampi**2*2D0 + p3p4*2D0) - 
     &             ampi**2*kp4*1D0*
     &              (m1**2 + ampi**2*2D0 + p3p4*2D0)*
     &              (-(1D0*(ame**2 - p1p2*1D0)*
     &                   (p2p3 - p2p4*1D0)) + 
     &                p1p3*(p2p3 + p2p4 + ame**2*2D0) - 
     &                p1p4*1D0*(p2p3 + p2p4 + ame**2*2D0)) + 
     &             p2p4*(2D0*
     &                 (ampi**2*
     &                    (ame**2*ampi**2 - 
     &                     1D0*(p1p3 - p1p4*1D0)*
     &                     (p2p3 - p2p4*1D0))*
     &                    (m1**2 + ampi**2*2D0) + 
     &                   ampi**2*p3p4*
     &                    (-((p1p3 - p1p4*1D0)*
     &                     (p2p3 - p2p4*1D0)*2D0) + 
     &                     ame**2*(m1**2 - ampi**2*2D0)) + 
     &                   ame**2*p3p4**2*
     &                    (-(m1**2*1D0) + ampi**2*4D0)) + 
     &                p1p2*2D0*
     &                 (ampi**4*(m1**2 + ampi**2*2D0) + 
     &                   p3p4*
     &                    (ampi**2*m1**2 - ampi**4*2D0 + 
     &                     p3p4*(-(m1**2*1D0) + ampi**2*4D0)))
     &                  + p3p4*x12*
     &                 (ampi**2*(ampi + m1)*(ampi - m1*1D0)*
     &                    2D0 + p3p4*(m1**2 - ampi**2*6D0))))
     &           + kp4*p2p3*
     &           (-(ampi**2*kp4*
     &                (ame**2 + p1p2 + p1p3 - p1p4*1D0)*2D0*
     &                (m1**2 + ampi**2*2D0 + p3p4*2D0)) + 
     &             ampi**2*
     &              (ame**2*ampi**2 - 
     &                1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0))
     &               *(m1**2 + ampi**2*2D0)*4D0 + 
     &             ame**2*p3p4**2*2D0*
     &              (-(m1**2*1D0) + ampi**2*4D0) + 
     &             ampi**2*p3p4*2D0*
     &              (ame**2*m1**2 - 
     &                (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*4D0)
     &               + p1p2*2D0*
     &              (ampi**2*m1**2*p3p4 + 
     &                ampi**4*2D0*(m1**2 + ampi**2*2D0) + 
     &                p3p4**2*(-(m1**2*1D0) + ampi**2*4D0)) + 
     &             p3p4*x12*
     &              (-(ampi**2*m1**2*3D0) + 
     &                p3p4*(m1**2 - ampi**2*8D0)))) + 
     &       kp2*(-(ampi**2*kp3**2*2D0*
     &             (p1p4*(ame**2 + p1p2 + p2p4 - p2p3*1D0) + 
     &               ame**2*kp4*2D0)*
     &             (m1**2 + ampi**2*2D0 + p3p4*2D0)) + 
     &          kp3*(ampi**2*kp4*2D0*
     &              (m1**2 + ampi**2*2D0 + p3p4*2D0)*
     &              (p1p4*
     &                 (p1p2 + p2p4 - ame**2*1D0 - p2p3*1D0)
     &                 + p1p3*
     &                 (ame**2 + p2p4 - p1p2*1D0 - p2p3*1D0)
     &                 + ame**2*(p2p4 - p2p3*1D0)*2D0) + 
     &             ame**2*ampi**2*kp4**2*
     &              (m1**2 + ampi**2*2D0 + p3p4*2D0)*4D0 + 
     &             p1p4*(-(ampi**2*kp2*
     &                   (ampi**2 + p1p4 - p1p3*1D0 - 
     &                     p3p4*1D0)*2D0*
     &                   (m1**2 + ampi**2*2D0 + p3p4*2D0)) + 
     &                ampi**2*
     &                 (ame**2*ampi**2 - 
     &                   1D0*(p1p3 - p1p4*1D0)*
     &                    (p2p3 - p2p4*1D0))*
     &                 (m1**2 + ampi**2*2D0)*4D0 - 
     &                p3p4*x12*1D0*(ampi**2 - p3p4*1D0)*
     &                 (-(m1**2*1D0) + ampi**2*4D0) + 
     &                ampi**2*p3p4*2D0*
     &                 (-((p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*
     &                     4D0) + 
     &                   ame**2*(-(m1**2*3D0) + ampi**2*4D0))
     &                 + ame**2*p3p4**2*2D0*
     &                 (m1**2 - ampi**2*8D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0)*2D0*
     &                 (ampi**2*2D0*(m1**2 + ampi**2*2D0) + 
     &                   p3p4*(-(m1**2*1D0) + ampi**2*8D0))))
     &           + kp4*p1p3*
     &           (-(ame**2*m1**2*p3p4**2*2D0) - 
     &             m1**2*p1p2*p3p4**2*2D0 + 
     &             ampi**2*kp4*
     &              (ame**2 + p1p2 + p2p3 - p2p4*1D0)*2D0*
     &              (m1**2 + ampi**2*2D0 + p3p4*2D0) + 
     &             ampi**2*kp2*
     &              (ampi**2 + p1p3 - p1p4*1D0 - p3p4*1D0)*
     &              2D0*(m1**2 + ampi**2*2D0 + p3p4*2D0) - 
     &             ampi**4*m1**2*p1p2*4D0 + 
     &             ampi**2*m1**2*p1p3*p2p3*4D0 - 
     &             ampi**2*m1**2*p1p4*p2p3*4D0 - 
     &             ampi**2*m1**2*p1p3*p2p4*4D0 + 
     &             ampi**2*m1**2*p1p4*p2p4*4D0 - 
     &             ame**2*ampi**4*(m1**2 + ampi**2*2D0)*4D0 + 
     &             p3p4*x12*(ampi**2 - p3p4*1D0)*
     &              (-(m1**2*1D0) + ampi**2*4D0) + 
     &             ame**2*ampi**2*m1**2*p3p4*6D0 + 
     &             ampi**2*m1**2*p1p2*p3p4*6D0 - 
     &             ampi**6*p1p2*8D0 + ampi**4*p1p3*p2p3*8D0 - 
     &             ampi**4*p1p4*p2p3*8D0 - 
     &             ampi**4*p1p3*p2p4*8D0 + 
     &             ampi**4*p1p4*p2p4*8D0 - 
     &             ame**2*ampi**4*p3p4*8D0 - 
     &             ampi**4*p1p2*p3p4*8D0 + 
     &             ampi**2*p1p3*p2p3*p3p4*8D0 - 
     &             ampi**2*p1p4*p2p3*p3p4*8D0 - 
     &             ampi**2*p1p3*p2p4*p3p4*8D0 + 
     &             ampi**2*p1p4*p2p4*p3p4*8D0 + 
     &             ame**2*ampi**2*p3p4**2*1.6D1 + 
     &             ampi**2*p1p2*p3p4**2*1.6D1))))/
     &   (ampi**2*kp1*kp2*kp3*kp4*(ame**2 + p1p2)*
     &     (ampi**4 - p3p4**2*1D0))
