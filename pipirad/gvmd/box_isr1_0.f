        elmat1 = -((C0aeex12m1em2*ep3*1D0*
     &       (kp1**2*(ampi**2*p1p2**2*(p2p3 - p2p4*1D0) + 
     &            ame**4*(ampi**2*p2p4 - p2p3*p3p4*1D0) + 
     &            ame**2*p1p3**2*(p2p3 - p2p4*1D0)*2D0 + 
     &            ame**2*p1p2*p2p3*(ampi**2 - p3p4*1D0)*2D0 + 
     &            p1p3*(ame**4*(p3p4 - ampi**2*1D0)*2D0 + 
     &               p1p2*
     &                (ame**2*(p3p4 - ampi**2*1D0) - 
     &                  p2p3**2*2D0 + p2p3*p2p4*2D0))) + 
     &         kp1*(ame**2*ampi**2*p1p2**2*p1p4 + 
     &            ampi**2*kp2*p1p2**2*p1p4 + 
     &            ampi**2*kp2*p1p2**2*p2p3 + 
     &            ame**4*p1p3*p1p4*p2p3 + 
     &            ame**4*p1p4*p2p3**2 + ame**6*ampi**2*p2p4 + 
     &            ame**4*ampi**2*kp2*p2p4 + 
     &            ame**4*p1p3*p2p3*p2p4 + 
     &            ame**2*p1p2*p1p3*p2p3*p2p4 + 
     &            ame**2*p1p2**2*p1p3*p3p4 + 
     &            ame**2*kp2*p1p2*p2p3*p3p4 + 
     &            ame**2*p1p2**2*p2p3*p3p4 - 
     &            ampi**2*kp2*p1p2**2*p1p3*1D0 - 
     &            ame**6*ampi**2*p1p4*1D0 - 
     &            ame**4*ampi**2*kp2*p1p4*1D0 - 
     &            ame**2*ampi**2*kp2*p1p2*p2p3*1D0 - 
     &            ame**2*p1p2*p1p4*p2p3**2*1D0 - 
     &            ame**2*ampi**2*p1p2**2*p2p4*1D0 - 
     &            ampi**2*kp2*p1p2**2*p2p4*1D0 - 
     &            ame**4*p1p3**2*p2p4*1D0 - 
     &            ame**2*p1p2*p1p3**2*p2p4*1D0 - 
     &            ame**6*p1p3*p3p4*1D0 - 
     &            ame**4*kp2*p1p3*p3p4*1D0 - 
     &            ame**6*p2p3*p3p4*1D0 + 
     &            ame**6*ampi**2*p1p3*2D0 + 
     &            ame**4*ampi**2*kp2*p1p3*2D0 - 
     &            ame**2*ampi**2*p1p2**2*p1p3*2D0 - 
     &            ame**4*p1p3**3*2D0 - 
     &            ame**2*kp2*p1p3**3*2D0 + 
     &            ame**4*p1p3**2*p1p4*2D0 + 
     &            ame**2*kp2*p1p3**2*p1p4*2D0 + 
     &            ame**2*kp2*p1p3**2*p2p3*2D0 + 
     &            kp2*p1p2*p1p3**2*p2p3*2D0 - 
     &            kp2*p1p2*p1p3*p1p4*p2p3*2D0 - 
     &            ame**4*p1p3*p2p3**2*2D0 + 
     &            ame**2*kp2*p2p3**3*2D0 - 
     &            ame**2*kp2*p1p3**2*p2p4*2D0 - 
     &            ame**2*kp2*p2p3**2*p2p4*2D0 - 
     &            kp4*1D0*(p1p2 + ame**2*2D0)*
     &             (-(ampi**2*p1p2**2*1D0) + 
     &               ame**2*
     &                (ame**2*ampi**2 - p1p3**2*1D0 - 
     &                  p2p3**2*1D0) + p1p2*p1p3*p2p3*2D0) + 
     &            kp3*(-(ampi**2*p1p2**3*1D0) + 
     &               p1p2**2*
     &                (-(ame**2*ampi**2*2D0) + p1p3*p2p3*2D0)
     &                + ame**4*
     &                (p1p4*p2p3 + ame**2*ampi**2*2D0 - 
     &                  p2p3*p2p4*2D0 + 
     &                  p1p3*(p2p4 - p1p4*2D0 - p2p3*2D0)) + 
     &               ame**2*p1p2*
     &                (ame**2*ampi**2 - p2p3*p2p4*1D0 + 
     &                  p1p4*p2p3*2D0 - 
     &                  p1p3*1D0*(p1p4 - p2p4*2D0))) + 
     &            ame**2*ampi**2*kp2*p1p2*p1p3*3D0 - 
     &            ame**2*p1p2*p1p3*p1p4*p2p3*3D0 - 
     &            ame**2*kp2*p1p2*p1p3*p3p4*3D0 + 
     &            ame**4*kp2*p2p3*p3p4*3D0 - 
     &            ame**4*ampi**2*kp2*p2p3*4D0 + 
     &            ame**2*p1p2*p1p3**2*p2p3*4D0 - 
     &            kp2*p1p2*p1p3*p2p3**2*4D0 + 
     &            kp2*p1p2*p1p3*p2p3*p2p4*4D0) + 
     &         kp2*(ame**2*ampi**2*p1p2**2*p1p4 + 
     &            ampi**2*kp2*p1p2**2*p1p4 + 
     &            ame**2*ampi**2*kp2*p1p2*p2p3 + 
     &            ame**4*p1p3*p1p4*p2p3 + 
     &            ame**4*p1p4*p2p3**2 + ame**6*ampi**2*p2p4 + 
     &            ame**2*p1p2*p1p3**2*p2p4 + 
     &            ame**6*p1p3*p3p4 + 
     &            ame**2*p1p2**2*p2p3*p3p4 - 
     &            ame**6*ampi**2*p1p4*1D0 - 
     &            ame**4*ampi**2*kp2*p1p4*1D0 - 
     &            ame**2*p1p2*p1p3*p1p4*p2p3*1D0 - 
     &            ame**2*p1p2*p1p4*p2p3**2*1D0 - 
     &            ame**2*ampi**2*p1p2**2*p2p4*1D0 - 
     &            ame**4*p1p3**2*p2p4*1D0 - 
     &            ame**4*p1p3*p2p3*p2p4*1D0 - 
     &            ame**2*p1p2*p1p3*p2p3*p2p4*1D0 - 
     &            ame**4*kp2*p1p3*p3p4*1D0 - 
     &            ame**2*p1p2**2*p1p3*p3p4*1D0 - 
     &            ame**6*p2p3*p3p4*1D0 - 
     &            ame**2*kp2*p1p2*p2p3*p3p4*1D0 + 
     &            ame**4*ampi**2*kp2*p1p3*2D0 - 
     &            ame**4*ampi**2*p1p2*p1p3*2D0 + 
     &            ampi**2*p1p2**3*p1p3*2D0 - 
     &            ame**2*kp2*p1p3**3*2D0 + 
     &            ame**2*p1p2*p1p3**3*2D0 + 
     &            ame**2*kp2*p1p3**2*p1p4*2D0 - 
     &            ame**2*p1p2*p1p3**2*p1p4*2D0 + 
     &            p1p2**2*p1p3*p1p4*p2p3*2D0 - 
     &            ame**2*kp2*p1p3*p2p3**2*2D0 + 
     &            ame**2*p1p2*p1p3*p2p3**2*2D0 + 
     &            ame**2*kp2*p1p4*p2p3**2*2D0 + 
     &            p1p2**2*p1p3**2*p2p4*2D0 + 
     &            ame**4*p1p2*p1p3*p3p4*2D0 + 
     &            kp2*p1p2**2*p1p3*p3p4*2D0 - 
     &            p1p2**3*p1p3*p3p4*2D0 + 
     &            kp4*(p1p2 + ame**2*2D0)*
     &             (-(ampi**2*p1p2**2*1D0) + 
     &               ame**2*
     &                (ame**2*ampi**2 - p1p3**2*1D0 - 
     &                  p2p3**2*1D0) + p1p2*p1p3*p2p3*2D0) - 
     &            ampi**2*kp2*p1p2**2*p1p3*3D0 + 
     &            kp2*p1p2*p1p3**2*p2p3*4D0 - 
     &            p1p2**2*p1p3**2*p2p3*4D0 - 
     &            kp2*p1p2*p1p3*p1p4*p2p3*4D0 + 
     &            kp3*(ampi**2*p1p2**3 + 
     &               p1p2**2*(ame**2*ampi**2 - p1p3*p2p4*1D0)*
     &                2D0 + 
     &               ame**4*
     &                (p1p3*p2p4 - p1p4*p2p3*1D0 - 
     &                  ame**2*ampi**2*2D0 + p1p3**2*2D0 + 
     &                  p2p3**2*2D0) + 
     &               ame**2*p1p2*
     &                (p2p3*p2p4 - ame**2*ampi**2*1D0 + 
     &                  p1p3*(p1p4 - p2p3*4D0))))))/
     &     (kp1*kp2*(ampi**2 + p3p4)*
     &       (ampi**2*p1p2**2 + 
     &         ame**2*(p1p3**2 + p2p3**2 - 
     &            ame**2*ampi**2*1D0) - p1p2*p1p3*p2p3*2D0)))
     &   - (C0aepx13em1p*ep3*1D0*
     &     (kp1**2*(ame**2*ampi**2*
     &           (p2p3**2 - p2p3*p2p4*1D0 + 
     &             ame**2*(p3p4 - ampi**2*1D0)) + 
     &          ame**2*p1p3*p2p3*(ampi**2 - p3p4*1D0)*2D0 + 
     &          p1p3**2*(ame**2*(ampi**2 - p3p4*1D0) - 
     &             p2p3**2*2D0 + p2p3*p2p4*2D0) + 
     &          ampi**2*p1p2*
     &           (p1p3*(p2p3 - p2p4*1D0) + 
     &             ame**2*(p3p4 - ampi**2*1D0)*2D0)) + 
     &       kp1*(ampi**2*p1p2**2*
     &           (kp4*p1p3 - kp3*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0)) + 
     &          p1p2*(ame**2*ampi**2*p1p4*p2p3 - 
     &             p1p3**2*1D0*(ame**2*p3p4 + kp4*p2p3*2D0) + 
     &             ame**2*p1p3*
     &              (ampi**2*p1p4 + p2p3*(p3p4 - ampi**2*2D0))
     &               + kp3*
     &              (ame**2*ampi**2*p1p3*2D0 + 
     &                p1p3**2*p2p3*2D0 + 
     &                ame**2*ampi**2*
     &                 (p2p4 - p2p3*1D0 - p1p4*2D0))) + 
     &          ame**2*(ame**2*ampi**2*p1p3**2 + 
     &             ame**2*ampi**2*p2p3**2 + p1p3**3*p2p4 + 
     &             p1p3**2*p2p3*p2p4 + ame**4*ampi**2*p3p4 + 
     &             ame**2*p1p3*p2p3*p3p4 - 
     &             ame**4*ampi**4*1D0 - 
     &             p1p3**2*p1p4*p2p3*1D0 - 
     &             p1p3*p1p4*p2p3**2*1D0 - 
     &             ame**2*ampi**2*p1p3*p2p4*1D0 - 
     &             ame**2*ampi**2*p2p3*p2p4*1D0 - 
     &             ame**2*p1p3**2*p3p4*1D0 + 
     &             kp4*p1p3*
     &              (p1p3**2 + p2p3**2 - ame**2*ampi**2*1D0)
     &              + kp3*
     &              (-(p1p3**3*2D0) + 
     &                p1p3*
     &                 (-(p2p3*p2p4*1D0) + 
     &                   ame**2*ampi**2*2D0 + p1p4*p2p3*2D0)
     &                 + p1p3**2*(p1p4 - p2p4*2D0) - 
     &                ame**2*ampi**2*1D0*
     &                 (p1p4 + p2p3*2D0 - p2p4*2D0))) + 
     &          kp2*(ampi**2*p1p2*
     &              (-(p1p3**2*1D0) + 
     &                p1p3*(p1p4 + p2p3 - p2p4*1D0) + 
     &                ame**2*(ampi**2 - p3p4*1D0)) + 
     &             p1p3**3*p2p3*2D0 + 
     &             ame**2*p1p3*p2p3*(p3p4 - ampi**2*2D0) + 
     &             p1p3**2*
     &              (-(p1p4*p2p3*2D0) + 
     &                ame**2*(p3p4 - ampi**2*1D0)*3D0) + 
     &             ame**2*ampi**2*
     &              (p1p4*p2p3 + p2p3*p2p4 - p2p3**2*1D0 + 
     &                ame**2*(ampi**2 - p3p4*1D0)*3D0))) + 
     &       kp2*(p1p2**2*
     &           (ampi**2*kp3*p1p3 - ampi**2*kp4*p1p3*1D0 + 
     &             (ampi**2 - p3p4*1D0)*
     &              (ame**2*ampi**2 - p1p3**2*2D0)) + 
     &          kp2*(ame**2*p2p3*
     &              (p1p3*p3p4 - ampi**2*p1p4*1D0) + 
     &             p1p2*(ampi**2*p1p3*p1p4 + 
     &                ame**2*ampi**2*(p3p4 - ampi**2*1D0) + 
     &                p1p3**2*(ampi**2 - p3p4*2D0))) - 
     &          ame**2*1D0*
     &           (ame**4*ampi**4 + p1p3**3*p2p4 + 
     &             ame**2*ampi**2*p2p3*p2p4 + 
     &             ame**2*p1p3**2*p3p4 + 
     &             ame**2*p1p3*p2p3*p3p4 - 
     &             p1p3**2*p1p4*p2p3*1D0 - 
     &             ame**2*ampi**2*p2p3**2*1D0 - 
     &             p1p3*p1p4*p2p3**2*1D0 - 
     &             ame**2*ampi**2*p1p3*p2p4*1D0 - 
     &             p1p3**2*p2p3*p2p4*1D0 - 
     &             ame**4*ampi**2*p3p4*1D0 + 
     &             kp4*p1p3*
     &              (p1p3**2 + p2p3**2 - ame**2*ampi**2*1D0)
     &              + p1p3**4*2D0 + 
     &             ame**2*ampi**2*p1p3*p1p4*2D0 - 
     &             p1p3**3*p1p4*2D0 + p1p3**2*p2p3**2*2D0 + 
     &             kp3*(p1p3**2*p1p4 - 
     &                ame**2*ampi**2*p1p4*1D0 - p1p3**3*2D0 + 
     &                p1p3*
     &                 (p2p3*p2p4 + ame**2*ampi**2*2D0 - 
     &                   p2p3**2*2D0)) - 
     &             ame**2*ampi**2*p1p3**2*3D0) + 
     &          p1p2*(ame**2*ampi**2*
     &              (p1p4*p2p3 + kp3*(p2p3 - p2p4*1D0)) - 
     &             ame**2*p1p3*1D0*
     &              (ampi**2*p1p4 - ampi**2*p2p4*2D0 + 
     &                p2p3*(p3p4 + ampi**2*2D0)) + 
     &             p1p3**3*(-(p2p4*2D0) + p2p3*4D0) + 
     &             p1p3**2*
     &              (ame**2*p3p4 + kp4*p2p3*2D0 - 
     &                p1p4*p2p3*2D0 + 
     &                kp3*(p2p4*2D0 - p2p3*4D0))))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*
     &     (ampi**2*p1p2**2 + 
     &       ame**2*(p1p3**2 + p2p3**2 - 
     &          ame**2*ampi**2*1D0) - p1p2*p1p3*p2p3*2D0)) - 
     &  (C0ax12px45m2m1p*ep3*1D0*
     &     (kp1**2*(ampi**2*p1p2**2*(p2p3 - p2p4*1D0) + 
     &          p1p3**2*(-(ame**2*1D0*
     &                (ampi**2 + p2p4 - p3p4*1D0)) + 
     &             p2p3**2*2D0 + p2p3*(ame**2 - p2p4*2D0)) + 
     &          ame**2*(p2p3**3 + 
     &             ampi**2*p2p3*(p2p4 - ame**2*1D0) - 
     &             p2p3**2*1D0*
     &              (p2p4 - p3p4*2D0 + ampi**2*3D0) + 
     &             ame**2*ampi**2*
     &              (p2p4 + ampi**2*3D0 - p3p4*3D0)) + 
     &          p1p2*(-(p1p3*1D0*(p2p3 - p2p4*1D0)*
     &                (ampi**2 + p2p3*2D0)) + 
     &             ampi**2*
     &              (p2p3*p2p4 - p2p3**2*1D0 + 
     &                ame**2*(ampi**2 - p3p4*1D0)*3D0)) + 
     &          p1p3*(ame**2*ampi**2*p2p4 + p2p3**3*2D0 - 
     &             p2p3**2*p2p4*2D0 + 
     &             ame**2*p2p3*(p3p4*3D0 - ampi**2*4D0))) + 
     &       kp1*(ampi**2*kp4*p1p2**3 + ame**6*ampi**2*p1p3 + 
     &          ame**4*ampi**2*kp4*p1p3 + 
     &          ame**2*kp4*p1p2*p1p3**2 + 
     &          ame**2*ampi**2*p1p2**2*p1p4 + 
     &          ame**4*p1p3**2*p1p4 + 
     &          ame**4*ampi**2*kp4*p2p3 + 
     &          ame**2*ampi**2*p1p2**2*p2p3 + 
     &          ame**4*p1p3**2*p2p3 + 
     &          ame**2*p1p3**2*p1p4*p2p3 + 
     &          ame**2*kp4*p1p2*p2p3**2 + 
     &          ame**4*p1p4*p2p3**2 + ame**4*p2p3**3 + 
     &          ame**2*p1p4*p2p3**3 + ame**6*ampi**2*p2p4 + 
     &          ame**4*ampi**2*p1p3*p2p4 + 
     &          ame**2*ampi**2*p1p2*p1p3*p2p4 + 
     &          ame**4*ampi**2*p2p3*p2p4 + 
     &          ame**2*ampi**2*p1p2*p2p3*p2p4 + 
     &          ame**4*p1p3**2*p3p4 + 
     &          ame**2*p1p2*p1p3**2*p3p4 - 
     &          ame**4*ampi**2*kp4*p1p2*1D0 - 
     &          ame**2*ampi**2*p1p2**2*p1p3*1D0 - 
     &          ampi**2*kp4*p1p2**2*p1p3*1D0 - 
     &          ame**4*p1p3**3*1D0 - ame**2*kp4*p1p3**3*1D0 - 
     &          ame**6*ampi**2*p1p4*1D0 - 
     &          ame**4*ampi**2*p1p3*p1p4*1D0 - 
     &          ame**2*ampi**2*p1p2*p1p3*p1p4*1D0 - 
     &          ame**6*ampi**2*p2p3*1D0 - 
     &          ampi**2*kp4*p1p2**2*p2p3*1D0 - 
     &          ame**2*kp4*p1p3**2*p2p3*1D0 - 
     &          ame**4*ampi**2*p1p4*p2p3*1D0 - 
     &          ame**2*ampi**2*p1p2*p1p4*p2p3*1D0 - 
     &          ame**4*p1p3*p2p3**2*1D0 - 
     &          ame**2*kp4*p1p3*p2p3**2*1D0 - 
     &          ame**2*kp4*p2p3**3*1D0 - 
     &          ame**2*ampi**2*p1p2**2*p2p4*1D0 - 
     &          ame**4*p1p3**2*p2p4*1D0 - 
     &          ame**2*p1p3**3*p2p4*1D0 - 
     &          ame**4*p2p3**2*p2p4*1D0 - 
     &          ame**2*p1p3*p2p3**2*p2p4*1D0 - 
     &          ame**4*p2p3**2*p3p4*1D0 - 
     &          ame**2*p1p2*p2p3**2*p3p4*1D0 - 
     &          ame**6*ampi**2*kp4*2D0 + 
     &          ame**2*ampi**2*kp4*p1p2**2*2D0 + 
     &          ame**4*kp4*p1p3**2*2D0 - 
     &          kp4*p1p2**2*p1p3*p2p3*2D0 + 
     &          ame**2*p1p2*p1p3**2*p2p3*2D0 + 
     &          kp4*p1p2*p1p3**2*p2p3*2D0 - 
     &          ame**2*p1p2*p1p3*p1p4*p2p3*2D0 + 
     &          ame**4*kp4*p2p3**2*2D0 - 
     &          ame**2*p1p2*p1p3*p2p3**2*2D0 + 
     &          kp4*p1p2*p1p3*p2p3**2*2D0 + 
     &          ame**2*p1p3*p1p4*p2p3**2*2D0 + 
     &          ame**2*p1p2*p1p3*p2p3*p2p4*2D0 - 
     &          ame**2*p1p3**2*p2p3*p2p4*2D0 - 
     &          ame**2*kp4*p1p2*p1p3*p2p3*4D0 + 
     &          kp3*(-(ampi**2*p1p2**3*1D0) + 
     &             p1p2**2*
     &              (ampi**2*(p2p3 - ame**2*2D0) + 
     &                p1p3*(ampi**2 + p2p3*2D0)) + 
     &             ame**2*
     &              (p2p3**2*p2p4 + ame**4*ampi**2*2D0 + 
     &                p1p3**3*2D0 + ame**2*ampi**2*p2p3*2D0 - 
     &                ame**2*p2p3**2*2D0 - 
     &                p1p3**2*1D0*
     &                 (p1p4 + ame**2*2D0 - p2p3*2D0 - 
     &                   p2p4*2D0) - 
     &                ame**2*ampi**2*p2p4*3D0 + 
     &                p1p4*
     &                 (-(p2p3**2*2D0) + ame**2*ampi**2*3D0)
     &                 + p1p3*
     &                 (-(p1p4*p2p3*3D0) + p2p3*p2p4*3D0 - 
     &                   ame**2*ampi**2*4D0)) + 
     &             p1p2*(-(p1p3**2*1D0*(ame**2 + p2p3*2D0)) + 
     &                ame**2*
     &                 (-(p2p3**2*1D0) + ampi**2*p1p4*3D0 + 
     &                   ampi**2*p2p3*3D0 + 
     &                   ampi**2*(ame**2 - p2p4*3D0)) + 
     &                p1p3*
     &                 (-(p2p3**2*2D0) - ame**2*ampi**2*3D0 + 
     &                   ame**2*p2p3*4D0))) - 
     &          kp2*1D0*(ampi**2*p1p2**2*
     &              (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0) + 
     &             p1p3**3*(ame**2 + p2p3*2D0) - 
     &             p1p3**2*1D0*
     &              (ame**2*p2p3 - p2p3**2*2D0 + 
     &                p1p4*(ame**2 + p2p3*2D0) + 
     &                ame**2*
     &                 (-(p2p4*1D0) - p3p4*3D0 + ampi**2*4D0))
     &               + ame**2*
     &              (ampi**2*p2p3*(ame**2 + p2p4) - 
     &                p2p3**3*1D0 + 
     &                p1p4*
     &                 (ame**2*ampi**2 + ampi**2*p2p3 - 
     &                   p2p3**2*1D0) + 
     &                p2p3**2*(p2p4 + p3p4 - ampi**2*2D0) + 
     &                ame**2*ampi**2*
     &                 (-(p2p4*1D0) + ampi**2*4D0 - p3p4*4D0))
     &               + p1p2*
     &              (-(p1p3**2*1D0*(ampi**2 + p2p3*2D0)) + 
     &                p1p3*
     &                 (-(ampi**2*p2p4*1D0) + p2p3**2*2D0 - 
     &                   p2p3*p2p4*2D0 + 
     &                   p1p4*(ampi**2 + p2p3*2D0)) + 
     &                ampi**2*
     &                 (p1p4*p2p3 + p2p3*p2p4 - p2p3**2*1D0 + 
     &                   ame**2*(ampi**2 - p3p4*1D0)*4D0)) + 
     &             p1p3*(p1p4*
     &                 (ame**2*ampi**2 - p2p3**2*2D0) + 
     &                ame**2*
     &                 (p2p3**2 - 
     &                   ampi**2*(ame**2 + p2p4)*1D0 - 
     &                   p2p3*(ampi**2 - p3p4*1D0)*4D0)))) + 
     &       kp2*(ame**4*ampi**2*kp4*p1p2 + 
     &          ame**6*ampi**2*p1p3 + 
     &          ampi**2*kp4*p1p2**2*p1p3 + 
     &          ame**2*kp4*p1p3**3 + 
     &          ame**2*ampi**2*p1p2**2*p1p4 + 
     &          ame**4*p1p3**2*p1p4 + 
     &          ame**2*ampi**2*p1p2**2*p2p3 + 
     &          ampi**2*kp4*p1p2**2*p2p3 + 
     &          ame**4*p1p3**2*p2p3 + 
     &          ame**2*kp4*p1p3**2*p2p3 + 
     &          ame**4*ampi**2*p1p4*p2p3 + 
     &          ame**2*kp4*p1p3*p2p3**2 + 
     &          ame**4*p1p4*p2p3**2 + ame**4*p2p3**3 + 
     &          ame**2*kp4*p2p3**3 + ame**6*ampi**2*p2p4 + 
     &          ame**2*p1p3**3*p2p4 + 
     &          ame**4*ampi**2*p2p3*p2p4 + 
     &          ame**2*ampi**2*p1p2*p2p3*p2p4 + 
     &          ame**4*p1p3**2*p3p4 + ame**4*p2p3**2*p3p4 + 
     &          ame**2*p1p2*p2p3**2*p3p4 - 
     &          ampi**2*kp4*p1p2**3*1D0 - 
     &          ame**4*ampi**2*kp4*p1p3*1D0 - 
     &          ame**2*ampi**2*p1p2**2*p1p3*1D0 - 
     &          ame**2*kp4*p1p2*p1p3**2*1D0 - 
     &          ame**4*p1p3**3*1D0 - 
     &          ame**6*ampi**2*p1p4*1D0 - 
     &          ame**6*ampi**2*p2p3*1D0 - 
     &          ame**4*ampi**2*kp4*p2p3*1D0 - 
     &          ame**2*ampi**2*p1p2*p1p4*p2p3*1D0 - 
     &          ame**2*kp4*p1p2*p2p3**2*1D0 - 
     &          ame**4*p1p3*p2p3**2*1D0 - 
     &          ame**2*p1p4*p2p3**3*1D0 - 
     &          ame**2*ampi**2*p1p2**2*p2p4*1D0 - 
     &          ame**4*ampi**2*p1p3*p2p4*1D0 - 
     &          ame**4*p1p3**2*p2p4*1D0 - 
     &          ame**4*p2p3**2*p2p4*1D0 - 
     &          ame**2*p1p3*p2p3**2*p2p4*1D0 - 
     &          ame**2*p1p2*p1p3**2*p3p4*1D0 + 
     &          ame**6*ampi**4*2D0 + ame**6*ampi**2*kp4*2D0 + 
     &          ame**4*ampi**4*p1p2*2D0 - 
     &          ame**2*ampi**4*p1p2**2*2D0 - 
     &          ame**2*ampi**2*kp4*p1p2**2*2D0 - 
     &          ampi**4*p1p2**3*2D0 - 
     &          ame**4*kp4*p1p3**2*2D0 - 
     &          ame**2*ampi**2*p1p2*p1p3**2*2D0 + 
     &          ampi**2*p1p2**2*p1p3**2*2D0 + 
     &          ame**2*p1p3**4*2D0 - 
     &          ame**2*p1p3**3*p1p4*2D0 - 
     &          ame**4*ampi**2*p1p3*p2p3*2D0 + 
     &          kp4*p1p2**2*p1p3*p2p3*2D0 + 
     &          ame**2*p1p2*p1p3**2*p2p3*2D0 - 
     &          kp4*p1p2*p1p3**2*p2p3*2D0 + 
     &          ame**2*p1p3**3*p2p3*2D0 - 
     &          ampi**2*p1p2**2*p1p4*p2p3*2D0 - 
     &          ame**2*p1p2*p1p3*p1p4*p2p3*2D0 + 
     &          p1p2*p1p3**2*p1p4*p2p3*2D0 - 
     &          ame**4*ampi**2*p2p3**2*2D0 - 
     &          ame**4*kp4*p2p3**2*2D0 - 
     &          ame**2*ampi**2*p1p2*p2p3**2*2D0 - 
     &          ame**2*p1p2*p1p3*p2p3**2*2D0 - 
     &          kp4*p1p2*p1p3*p2p3**2*2D0 + 
     &          ame**2*p1p3**2*p2p3**2*2D0 - 
     &          ame**2*p1p3*p1p4*p2p3**2*2D0 + 
     &          p1p2*p1p3*p1p4*p2p3**2*2D0 + 
     &          ame**2*p1p3*p2p3**3*2D0 - 
     &          ampi**2*p1p2**2*p1p3*p2p4*2D0 + 
     &          p1p2*p1p3**3*p2p4*2D0 + 
     &          ame**2*p1p2*p1p3*p2p3*p2p4*2D0 + 
     &          p1p2*p1p3**2*p2p3*p2p4*2D0 - 
     &          ame**6*ampi**2*p3p4*2D0 - 
     &          ame**4*ampi**2*p1p2*p3p4*2D0 + 
     &          ame**2*ampi**2*p1p2**2*p3p4*2D0 + 
     &          ampi**2*p1p2**3*p3p4*2D0 - 
     &          p1p2**2*p1p3**2*p3p4*2D0 + 
     &          ame**4*p1p3*p2p3*p3p4*2D0 - 
     &          p1p2**2*p1p3*p2p3*p3p4*2D0 + 
     &          ame**4*ampi**2*p1p3*p1p4*3D0 + 
     &          ame**2*ampi**2*p1p2*p1p3*p1p4*3D0 - 
     &          ame**2*p1p3**2*p1p4*p2p3*3D0 - 
     &          ame**2*ampi**2*p1p2*p1p3*p2p4*3D0 + 
     &          kp2*(ame**2*
     &              (p1p3**2*(ampi**2 + p1p4) - p1p3**3*1D0 + 
     &                p1p4*
     &                 (ampi**2*p2p3 + p2p3**2 - 
     &                   ame**2*ampi**2*1D0) - 
     &                1D0*(ame**2*ampi**2 - p2p3**2*1D0)*
     &                 (ampi**2 - p3p4*1D0) + 
     &                p1p3*
     &                 (ame**2*ampi**2 - ampi**2*p1p4*1D0 - 
     &                   p2p3**2*1D0 - p2p3*p3p4*1D0)) + 
     &             ampi**2*p1p2**2*
     &              (p1p4 - p1p3*1D0 + ampi**2*2D0 - p3p4*2D0)
     &               + p1p2*
     &              (ampi**2*
     &                 (p1p4*p2p3 + 
     &                   ame**2*(ampi**2 - p3p4*1D0)) + 
     &                p1p3**2*
     &                 (-(ampi**2*1D0) + p2p3*2D0 + p3p4*2D0)
     &                 - p1p3*1D0*
     &                 (p1p4*(ampi**2 + p2p3*2D0) + 
     &                   p2p3*(-(p3p4*2D0) + ampi**2*3D0))))
     &           - ame**4*ampi**2*p1p3**2*4D0 + 
     &          ame**2*ampi**2*p1p2*p1p3*p2p3*4D0 + 
     &          ame**2*kp4*p1p2*p1p3*p2p3*4D0 - 
     &          p1p2*p1p3**3*p2p3*4D0 - 
     &          p1p2*p1p3**2*p2p3**2*4D0 + 
     &          kp3*(ampi**2*p1p2**3 + 
     &             ame**2*
     &              (-(ame**2*ampi**2*p1p4*1D0) - 
     &                p1p3**3*2D0 + 
     &                p1p3**2*
     &                 (p1p4 + ame**2*2D0 - p2p3*2D0) - 
     &                1D0*(ame**2*ampi**2 - p2p3**2*1D0)*
     &                 (p2p4 + ame**2*2D0 - p2p3*2D0) + 
     &                p1p3*
     &                 (p1p4*p2p3 + p2p3*p2p4 + 
     &                   ame**2*ampi**2*2D0 - p2p3**2*2D0)) + 
     &             p1p2**2*
     &              (-(p1p3*1D0*(ampi**2 + p2p3*2D0)) + 
     &                ampi**2*((ame**2 + p2p4)*2D0 - p2p3*3D0)
     &                ) + 
     &             p1p2*(ame**2*
     &                 (p2p3**2 - ampi**2*p1p4*1D0 - 
     &                   ampi**2*p2p3*1D0 + 
     &                   ampi**2*(p2p4 - ame**2*1D0)) + 
     &                p1p3**2*
     &                 (ame**2 - p2p4*2D0 + p2p3*4D0) + 
     &                p1p3*
     &                 (ame**2*ampi**2 - 
     &                   p2p3*2D0*(p2p4 + ame**2*2D0) + 
     &                   p2p3**2*4D0))) + 
     &          ampi**2*p1p2**2*p1p3*p2p3*6D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*
     &     (ampi**2*p1p2**2 + 
     &       ame**2*(p1p3**2 + p2p3**2 - 
     &          ame**2*ampi**2*1D0) - p1p2*p1p3*p2p3*2D0)) + 
     &  (C0aex13x45m2ep*ep3*
     &     (kp1**2*(-(ame**4*p2p3*p3p4*1D0) + 
     &          ampi**2*p1p2**2*(p2p3 - p2p4*1D0) + 
     &          ame**2*p2p3**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &          ame**2*p1p3**2*(p2p3 - p2p4*1D0)*2D0 + 
     &          ame**4*ampi**2*
     &           (p2p4 + ampi**2*2D0 - p3p4*2D0) + 
     &          p1p3*(ame**2*ampi**2*p2p4 + p2p3**3*2D0 - 
     &             p2p3**2*p2p4*2D0 + 
     &             ame**4*(p3p4 - ampi**2*1D0)*2D0 + 
     &             ame**2*p2p3*(p3p4 - ampi**2*2D0)) + 
     &          p1p2*(-(ampi**2*p2p3**2*1D0) + 
     &             ame**2*ampi**2*(ampi**2 - p3p4*1D0) + 
     &             p1p3*(ame**2*(p3p4 - ampi**2*1D0) - 
     &                p2p3**2*2D0 + p2p3*p2p4*2D0) + 
     &             p2p3*(ampi**2*p2p4 + 
     &                ame**2*(ampi**2 - p3p4*1D0)*2D0))) + 
     &       kp1*(ame**2*ampi**4*p1p2**2 + 
     &          ampi**2*kp4*p1p2**3 + 
     &          ame**4*ampi**2*p1p3**2 + 
     &          ame**2*kp4*p1p2*p1p3**2 + 
     &          ame**2*ampi**2*p1p2**2*p1p4 + 
     &          ame**4*ampi**2*kp4*p2p3 + 
     &          ame**4*p1p3*p1p4*p2p3 + 
     &          ame**4*ampi**2*p2p3**2 + 
     &          ame**2*kp4*p1p2*p2p3**2 + 
     &          ame**4*p1p4*p2p3**2 + 
     &          ame**2*p1p3*p1p4*p2p3**2 + 
     &          ame**2*p1p4*p2p3**3 + ame**6*ampi**2*p2p4 + 
     &          ame**2*ampi**2*p1p2*p1p3*p2p4 + 
     &          ame**2*ampi**2*p1p2*p2p3*p2p4 + 
     &          ame**4*p1p3*p2p3*p2p4 + 
     &          ame**2*p1p2*p1p3*p2p3*p2p4 + 
     &          ame**6*ampi**2*p3p4 + 
     &          ame**2*p1p2**2*p1p3*p3p4 + 
     &          ame**2*p1p2**2*p2p3*p3p4 + 
     &          ame**4*p1p3*p2p3*p3p4 + 
     &          ame**2*p1p2*p1p3*p2p3*p3p4 - 
     &          ame**6*ampi**4*1D0 - 
     &          ame**4*ampi**2*kp4*p1p2*1D0 - 
     &          ame**6*ampi**2*p1p4*1D0 - 
     &          ame**4*ampi**2*p1p3*p1p4*1D0 - 
     &          ampi**2*kp4*p1p2**2*p2p3*1D0 - 
     &          ame**2*kp4*p1p3**2*p2p3*1D0 - 
     &          ame**4*ampi**2*p1p4*p2p3*1D0 - 
     &          ame**2*p1p2*p1p4*p2p3**2*1D0 - 
     &          ame**2*kp4*p2p3**3*1D0 - 
     &          ame**2*ampi**2*p1p2**2*p2p4*1D0 - 
     &          ame**4*p1p3**2*p2p4*1D0 - 
     &          ame**2*p1p2*p1p3**2*p2p4*1D0 - 
     &          ame**2*p1p3**2*p2p3*p2p4*1D0 - 
     &          ame**2*p1p3*p2p3**2*p2p4*1D0 - 
     &          ame**2*ampi**2*p1p2**2*p3p4*1D0 - 
     &          ame**6*p1p3*p3p4*1D0 - ame**6*p2p3*p3p4*1D0 - 
     &          ame**4*p2p3**2*p3p4*1D0 - 
     &          ame**2*p1p2*p2p3**2*p3p4*1D0 - 
     &          ame**6*ampi**2*kp4*2D0 + 
     &          ame**2*ampi**2*kp4*p1p2**2*2D0 + 
     &          ame**6*ampi**2*p1p3*2D0 - 
     &          ame**2*ampi**2*p1p2**2*p1p3*2D0 + 
     &          ame**4*kp4*p1p3**2*2D0 - ame**4*p1p3**3*2D0 + 
     &          ame**4*p1p3**2*p1p4*2D0 - 
     &          ame**2*ampi**2*p1p2*p1p3*p2p3*2D0 - 
     &          kp4*p1p2**2*p1p3*p2p3*2D0 + 
     &          ame**4*kp4*p2p3**2*2D0 - 
     &          ame**4*p1p3*p2p3**2*2D0 + 
     &          kp4*p1p2*p1p3*p2p3**2*2D0 + 
     &          kp3*(-(ampi**2*p1p2**3*1D0) + 
     &             p1p2**2*
     &              (-(ame**2*ampi**2*2D0) + 
     &                p2p3*(ampi**2 + p1p3*2D0)) + 
     &             ame**2*
     &              (p2p3**2*p2p4 - ame**2*ampi**2*p2p4*1D0 + 
     &                ame**4*ampi**2*2D0 + p1p3**2*p2p3*2D0 - 
     &                ame**2*p2p3*p2p4*2D0 + 
     &                p1p4*
     &                 (ame**2*p2p3 + ame**2*ampi**2*2D0 - 
     &                   p2p3**2*2D0) - 
     &                p1p3*1D0*
     &                 (p2p3*(ame**2 - p2p4*1D0)*2D0 + 
     &                   p1p4*(p2p3 + ame**2*2D0) + 
     &                   ame**2*(-(p2p4*1D0) + ampi**2*2D0)))
     &              + p1p2*
     &              (ame**2*
     &                 (p2p3*(-(p2p4*1D0) + ampi**2*2D0) + 
     &                   p1p4*(ampi**2 + p2p3*2D0) + 
     &                   ampi**2*(ame**2 - p2p4*2D0)) - 
     &                p1p3*1D0*
     &                 (ame**2*p1p4 + p2p3**2*2D0 + 
     &                   ame**2*(ampi**2 - p2p4*2D0)))) - 
     &          ame**2*p1p2*p1p3*p1p4*p2p3*3D0 - 
     &          ame**2*kp4*p1p2*p1p3*p2p3*4D0 + 
     &          ame**2*p1p2*p1p3**2*p2p3*4D0 + 
     &          kp2*(ampi**2*p1p2**2*
     &              (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0) - 
     &             ame**2*p1p3**3*2D0 + 
     &             p1p3**2*
     &              (ame**2*p1p4*2D0 + ame**2*p2p3*2D0 - 
     &                p2p3**2*2D0 + 
     &                ame**2*(ampi**2 - p2p4*2D0)) + 
     &             p1p3*(p1p4*
     &                 (-(ame**2*ampi**2*1D0) + p2p3**2*2D0)
     &                 + ame**2*
     &                 (ampi**2*p2p4 + 
     &                   ame**2*(-(p3p4*1D0) + ampi**2*2D0) + 
     &                   p2p3*(ampi**2*2D0 - p3p4*3D0))) + 
     &             ame**2*
     &              (-(ame**2*ampi**2*p1p4*1D0) + 
     &                ame**2*ampi**2*
     &                 (p2p4 + p3p4 - ampi**2*1D0) + 
     &                p2p3**3*2D0 + 
     &                p2p3**2*
     &                 (ampi**2 - p3p4*1D0 - p2p4*2D0) + 
     &                ame**2*p2p3*(p3p4*3D0 - ampi**2*4D0)) + 
     &             p1p2*(ampi**2*p2p3**2 - 
     &                p2p3*1D0*
     &                 (ampi**2*p1p4 + ampi**2*p2p4 + 
     &                   ame**2*(ampi**2 - p3p4*1D0)) + 
     &                p1p3**2*p2p3*2D0 + 
     &                ame**2*ampi**2*(p3p4 - ampi**2*1D0)*
     &                 3D0 + 
     &                p1p3*
     &                 (ame**2*(ampi**2 - p3p4*1D0)*3D0 - 
     &                   p2p3**2*4D0 + 
     &                   p2p3*(ampi**2 - p1p4*2D0 + p2p4*4D0))
     &                ))) + 
     &       kp2*(ame**6*ampi**4 + ame**4*ampi**2*kp4*p1p2 + 
     &          ame**2*ampi**2*p1p2**2*p1p4 + 
     &          ame**4*ampi**2*p1p3*p1p4 + 
     &          ampi**2*kp4*p1p2**2*p2p3 + 
     &          ame**2*kp4*p1p3**2*p2p3 + 
     &          ame**4*ampi**2*p1p4*p2p3 + 
     &          ame**4*p1p3*p1p4*p2p3 + ame**4*p1p4*p2p3**2 + 
     &          ame**2*kp4*p2p3**3 + ame**6*ampi**2*p2p4 + 
     &          ame**2*p1p2*p1p3**2*p2p4 + 
     &          ame**2*ampi**2*p1p2*p2p3*p2p4 + 
     &          ame**2*p1p3**2*p2p3*p2p4 + 
     &          ame**2*ampi**2*p1p2**2*p3p4 + 
     &          ame**6*p1p3*p3p4 + ame**2*p1p2**2*p2p3*p3p4 + 
     &          ame**4*p1p3*p2p3*p3p4 + ame**4*p2p3**2*p3p4 + 
     &          ame**2*p1p2*p2p3**2*p3p4 - 
     &          ame**2*ampi**4*p1p2**2*1D0 - 
     &          ampi**2*kp4*p1p2**3*1D0 - 
     &          ame**4*ampi**2*p1p3**2*1D0 - 
     &          ame**2*kp4*p1p2*p1p3**2*1D0 - 
     &          ame**6*ampi**2*p1p4*1D0 - 
     &          ame**4*ampi**2*kp4*p2p3*1D0 - 
     &          ame**2*p1p2*p1p3*p1p4*p2p3*1D0 - 
     &          ame**4*ampi**2*p2p3**2*1D0 - 
     &          ame**2*kp4*p1p2*p2p3**2*1D0 - 
     &          ame**2*p1p2*p1p4*p2p3**2*1D0 - 
     &          ame**2*p1p3*p1p4*p2p3**2*1D0 - 
     &          ame**2*p1p4*p2p3**3*1D0 - 
     &          ame**2*ampi**2*p1p2**2*p2p4*1D0 - 
     &          ame**2*ampi**2*p1p2*p1p3*p2p4*1D0 - 
     &          ame**4*p1p3**2*p2p4*1D0 - 
     &          ame**4*p1p3*p2p3*p2p4*1D0 - 
     &          ame**2*p1p2*p1p3*p2p3*p2p4*1D0 - 
     &          ame**2*p1p3*p2p3**2*p2p4*1D0 - 
     &          ame**6*ampi**2*p3p4*1D0 - 
     &          ame**2*p1p2**2*p1p3*p3p4*1D0 - 
     &          ame**6*p2p3*p3p4*1D0 - 
     &          ame**2*p1p2*p1p3*p2p3*p3p4*1D0 + 
     &          ame**6*ampi**2*kp4*2D0 + 
     &          ame**4*ampi**4*p1p2*2D0 - 
     &          ame**2*ampi**2*kp4*p1p2**2*2D0 - 
     &          ampi**4*p1p2**3*2D0 - 
     &          ame**4*ampi**2*p1p2*p1p3*2D0 + 
     &          ampi**2*p1p2**3*p1p3*2D0 - 
     &          ame**4*kp4*p1p3**2*2D0 - 
     &          ame**2*ampi**2*p1p2*p1p3**2*2D0 + 
     &          ame**2*p1p2*p1p3**3*2D0 + 
     &          ame**2*ampi**2*p1p2*p1p3*p1p4*2D0 - 
     &          ame**2*p1p2*p1p3**2*p1p4*2D0 - 
     &          ame**4*ampi**2*p1p3*p2p3*2D0 + 
     &          ame**2*ampi**2*p1p2*p1p3*p2p3*2D0 + 
     &          kp4*p1p2**2*p1p3*p2p3*2D0 + 
     &          ame**2*p1p3**3*p2p3*2D0 - 
     &          ampi**2*p1p2**2*p1p4*p2p3*2D0 + 
     &          p1p2**2*p1p3*p1p4*p2p3*2D0 - 
     &          ame**2*p1p3**2*p1p4*p2p3*2D0 - 
     &          ame**4*kp4*p2p3**2*2D0 - 
     &          ame**2*ampi**2*p1p2*p2p3**2*2D0 + 
     &          ame**2*p1p2*p1p3*p2p3**2*2D0 - 
     &          kp4*p1p2*p1p3*p2p3**2*2D0 + 
     &          p1p2*p1p3*p1p4*p2p3**2*2D0 + 
     &          ame**2*p1p3*p2p3**3*2D0 - 
     &          ampi**2*p1p2**2*p1p3*p2p4*2D0 + 
     &          p1p2**2*p1p3**2*p2p4*2D0 + 
     &          p1p2*p1p3**2*p2p3*p2p4*2D0 - 
     &          ame**4*ampi**2*p1p2*p3p4*2D0 + 
     &          ampi**2*p1p2**3*p3p4*2D0 + 
     &          ame**4*p1p2*p1p3*p3p4*2D0 - 
     &          p1p2**3*p1p3*p3p4*2D0 - 
     &          p1p2**2*p1p3*p2p3*p3p4*2D0 + 
     &          ame**2*kp4*p1p2*p1p3*p2p3*4D0 - 
     &          p1p2**2*p1p3**2*p2p3*4D0 - 
     &          p1p2*p1p3**2*p2p3**2*4D0 + 
     &          kp2*(-(ame**2*1D0*
     &                ((ame**2*ampi**2 - p2p3**2*1D0)*
     &                   (ampi**2 - p3p4*1D0) + p1p3**3*2D0 - 
     &                  p1p3**2*1D0*(ampi**2 + p1p4*2D0) + 
     &                  p1p4*(ame**2*ampi**2 - p2p3**2*2D0) + 
     &                  p1p3*
     &                   (ampi**2*p1p4 + p2p3**2*2D0 + 
     &                     ame**2*(p3p4 - ampi**2*2D0)))) + 
     &             p1p2**2*
     &              (ampi**2*
     &                 (p1p4 + ampi**2*2D0 - p3p4*2D0) + 
     &                p1p3*(p3p4*2D0 - ampi**2*3D0)) + 
     &             p1p2*p2p3*
     &              (ampi**2*p1p4 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p3**2*4D0 + 
     &                p1p3*(p3p4*2D0 - ampi**2*3D0 - p1p4*4D0)
     &                )) + 
     &          kp3*(ampi**2*p1p2**3 + 
     &             ame**2*
     &              (p1p3*(p1p4*p2p3 + ame**2*p2p4) - 
     &                p2p3**3*2D0 + 
     &                p1p3**2*(ame**2 - p2p3*1D0)*2D0 + 
     &                p2p3**2*(p2p4 + ame**2*2D0) - 
     &                ame**2*ampi**2*1D0*
     &                 (p2p4 + ame**2*2D0) + 
     &                ame**2*p2p3*(-(p1p4*1D0) + ampi**2*2D0))
     &               + p1p2**2*
     &              (ame**2*ampi**2*2D0 + 
     &                p2p4*(ampi**2 - p1p3*1D0)*2D0 - 
     &                ampi**2*p2p3*3D0) + 
     &             p1p2*(-(ame**2*1D0*
     &                   (ame**2*ampi**2 + ampi**2*p1p4 - 
     &                     p2p3*p2p4*1D0)) + 
     &                p1p3*
     &                 (ame**2*ampi**2 + ame**2*p1p4 - 
     &                   p2p3*2D0*(p2p4 + ame**2*2D0) + 
     &                   p2p3**2*4D0))) + 
     &          ampi**2*p1p2**2*p1p3*p2p3*6D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*
     &     (ampi**2*p1p2**2 + 
     &       ame**2*(p1p3**2 + p2p3**2 - 
     &          ame**2*ampi**2*1D0) - p1p2*p1p3*p2p3*2D0)) - 
     &  (D0aex12x45x13epem1m2p*ep3*1D0*
     &     (kp1**2*(ame**2*
     &           (p2p3*(ame**2*m12*p3p4 + 
     &                ampi**2*p2p4*(-(m22*1D0) + ame**2*2D0))
     &              + ame**2*ampi**2*
     &              (-(m12*p2p4*1D0) + 
     &                (ampi**2 - p3p4*1D0)*
     &                 (-(m22*1D0) + ame**2*2D0 - m12*2D0)) + 
     &             p2p3**2*
     &              (-(m12*p3p4*2D0) + 
     &                ampi**2*(m22 - ame**2*2D0 + m12*2D0)))
     &           + ame**2*p1p3**3*(p2p3 - p2p4*1D0)*4D0 + 
     &          ampi**2*p1p2**2*
     &           (p2p4*(m12 - p1p3*2D0) - 
     &             p2p3*1D0*(m12 - p1p3*2D0) + 
     &             ame**2*(ampi**2 - p3p4*1D0)*4D0) + 
     &          p1p3**2*(-(p2p3*2D0*
     &                (ame**2*m12 + 
     &                  p2p4*(-(m22*1D0) + ame**2*2D0))) + 
     &             ame**2*
     &              (m12*p2p4*2D0 - 
     &                1D0*(ampi**2 - p3p4*1D0)*
     &                 (-(m22*1D0) + ame**2*2D0)) + 
     &             p2p3**2*(-(m22*2D0) + ame**2*4D0)) + 
     &          p1p3*(ame**4*m12*(ampi**2 - p3p4*1D0)*2D0 + 
     &             p2p3**2*p2p4*2D0*(m12 - ame**2*2D0) + 
     &             ame**2*ampi**2*p2p4*
     &              (-(m12*1D0) + ame**2*4D0) + 
     &             p2p3**3*(-(m12*2D0) + ame**2*4D0) + 
     &             ame**2*p2p3*
     &              (ampi**2*2D0*(m12 + m22 - ame**2*4D0) + 
     &                p3p4*(-(m12*1D0) - m22*2D0 + ame**2*4D0)
     &                )) + 
     &          p1p2*(ampi**2*p2p3**2*(m12 - ame**2*2D0) + 
     &             p2p3*(ame**2*m12*(p3p4 - ampi**2*1D0)*
     &                 2D0 + 
     &                ampi**2*p2p4*(-(m12*1D0) + ame**2*2D0))
     &              + p1p3**2*
     &              (ame**2*(p3p4 - ampi**2*1D0)*2D0 - 
     &                p2p3**2*4D0 + p2p3*p2p4*4D0) + 
     &             ame**2*ampi**2*(ampi**2 - p3p4*1D0)*
     &              (-(m12*1D0) - m22*2D0 + ame**2*6D0) + 
     &             p1p3*(ame**2*m12*(ampi**2 - p3p4*1D0) + 
     &                m12*p2p3**2*2D0 + 
     &                ampi**2*p2p4*
     &                 (-(m22*1D0) + ame**2*2D0) + 
     &                p2p3*
     &                 (-(m12*p2p4*2D0) + ame**2*p3p4*4D0 + 
     &                   ampi**2*(m22 - ame**2*6D0))))) + 
     &       kp1*(ame**6*ampi**4*m12 + 
     &          ame**4*ampi**2*kp4*m12*p1p2 + 
     &          ame**2*ampi**4*m22*p1p2**2 + 
     &          ampi**2*kp4*m22*p1p2**2*p1p3 + 
     &          ame**4*ampi**2*m22*p1p3**2 + 
     &          ame**2*kp4*m22*p1p3**3 + 
     &          ame**6*ampi**2*m12*p1p4 + 
     &          ame**4*ampi**2*m12*p1p3*p1p4 + 
     &          ame**2*ampi**2*m22*p1p2*p1p3*p1p4 + 
     &          ampi**2*kp4*m12*p1p2**2*p2p3 + 
     &          ame**2*kp4*m12*p1p3**2*p2p3 + 
     &          ame**4*ampi**2*m12*p1p4*p2p3 + 
     &          ame**2*ampi**2*m22*p1p2*p1p4*p2p3 + 
     &          ame**4*ampi**2*m22*p2p3**2 + 
     &          ame**2*kp4*m22*p1p3*p2p3**2 + 
     &          ame**2*m12*p1p2*p1p4*p2p3**2 + 
     &          ame**2*kp4*m12*p2p3**3 + 
     &          ame**2*ampi**2*m12*p1p2**2*p2p4 + 
     &          ame**4*m12*p1p3**2*p2p4 + 
     &          ame**2*m12*p1p2*p1p3**2*p2p4 + 
     &          ame**2*m22*p1p3**3*p2p4 + 
     &          ame**2*m12*p1p3**2*p2p3*p2p4 + 
     &          ame**2*m22*p1p3**2*p2p3*p2p4 + 
     &          ame**2*m12*p1p3*p2p3**2*p2p4 + 
     &          ame**6*ampi**2*m22*p3p4 + 
     &          ame**2*ampi**2*m12*p1p2**2*p3p4 + 
     &          ame**6*m12*p1p3*p3p4 + ame**6*m12*p2p3*p3p4 + 
     &          ame**4*m22*p1p3*p2p3*p3p4 + 
     &          ame**2*m22*p1p2*p1p3*p2p3*p3p4 + 
     &          ame**4*m12*p2p3**2*p3p4 + 
     &          ame**2*m12*p1p2*p2p3**2*p3p4 - 
     &          ame**6*ampi**4*m22*1D0 - 
     &          ame**2*ampi**4*m12*p1p2**2*1D0 - 
     &          ampi**2*kp4*m12*p1p2**3*1D0 - 
     &          ame**4*ampi**2*kp4*m22*p1p3*1D0 - 
     &          ame**4*ampi**2*m12*p1p3**2*1D0 - 
     &          ame**2*kp4*m12*p1p2*p1p3**2*1D0 - 
     &          ame**2*ampi**2*m12*p1p2**2*p1p4*1D0 - 
     &          ame**4*ampi**2*kp4*m12*p2p3*1D0 - 
     &          ame**4*m12*p1p3*p1p4*p2p3*1D0 - 
     &          ame**2*m22*p1p3**2*p1p4*p2p3*1D0 - 
     &          ame**4*ampi**2*m12*p2p3**2*1D0 - 
     &          ame**2*kp4*m12*p1p2*p2p3**2*1D0 - 
     &          ame**4*m12*p1p4*p2p3**2*1D0 - 
     &          ame**2*m12*p1p3*p1p4*p2p3**2*1D0 - 
     &          ame**2*m22*p1p3*p1p4*p2p3**2*1D0 - 
     &          ame**2*m12*p1p4*p2p3**3*1D0 - 
     &          ame**6*ampi**2*m12*p2p4*1D0 - 
     &          ame**4*ampi**2*m22*p1p3*p2p4*1D0 - 
     &          ame**2*ampi**2*m12*p1p2*p1p3*p2p4*1D0 - 
     &          ame**4*ampi**2*m22*p2p3*p2p4*1D0 - 
     &          ame**2*ampi**2*m12*p1p2*p2p3*p2p4*1D0 - 
     &          ame**4*m12*p1p3*p2p3*p2p4*1D0 - 
     &          ame**2*m12*p1p2*p1p3*p2p3*p2p4*1D0 - 
     &          ame**6*ampi**2*m12*p3p4*1D0 - 
     &          ame**2*ampi**2*m22*p1p2**2*p3p4*1D0 - 
     &          ame**2*m12*p1p2**2*p1p3*p3p4*1D0 - 
     &          ame**4*m22*p1p3**2*p3p4*1D0 - 
     &          ame**2*m22*p1p2*p1p3**2*p3p4*1D0 - 
     &          ame**2*m12*p1p2**2*p2p3*p3p4*1D0 - 
     &          ame**4*m12*p1p3*p2p3*p3p4*1D0 - 
     &          ame**2*m12*p1p2*p1p3*p2p3*p3p4*1D0 + 
     &          ame**8*ampi**4*2D0 + 
     &          ame**6*ampi**2*kp4*m12*2D0 + 
     &          ame**6*ampi**4*p1p2*2D0 - 
     &          ame**4*ampi**4*p1p2**2*2D0 - 
     &          ame**2*ampi**2*kp4*m12*p1p2**2*2D0 - 
     &          ame**2*ampi**4*p1p2**3*2D0 - 
     &          ame**6*ampi**2*m12*p1p3*2D0 - 
     &          ame**4*ampi**2*kp4*p1p2*p1p3*2D0 + 
     &          ame**2*ampi**2*m12*p1p2**2*p1p3*2D0 + 
     &          ampi**2*kp4*p1p2**3*p1p3*2D0 + 
     &          ame**6*ampi**2*p1p3**2*2D0 - 
     &          ame**4*kp4*m12*p1p3**2*2D0 - 
     &          ame**4*ampi**2*p1p2*p1p3**2*2D0 + 
     &          ame**4*m12*p1p3**3*2D0 + 
     &          ame**2*kp4*p1p2*p1p3**3*2D0 - 
     &          ame**4*ampi**2*p1p2*p1p3*p1p4*2D0 + 
     &          ame**2*ampi**2*p1p2**2*p1p3*p1p4*2D0 - 
     &          ame**4*m12*p1p3**2*p1p4*2D0 + 
     &          ame**2*ampi**2*m12*p1p2*p1p3*p2p3*2D0 - 
     &          ame**2*ampi**2*m22*p1p2*p1p3*p2p3*2D0 + 
     &          kp4*m12*p1p2**2*p1p3*p2p3*2D0 - 
     &          kp4*m22*p1p2*p1p3**2*p2p3*2D0 - 
     &          ame**4*ampi**2*p1p2*p1p4*p2p3*2D0 - 
     &          ame**2*ampi**2*p1p2**2*p1p4*p2p3*2D0 + 
     &          ame**4*p1p3**2*p1p4*p2p3*2D0 - 
     &          ame**6*ampi**2*p2p3**2*2D0 - 
     &          ame**4*kp4*m12*p2p3**2*2D0 - 
     &          ame**4*ampi**2*p1p2*p2p3**2*2D0 + 
     &          ame**4*m12*p1p3*p2p3**2*2D0 + 
     &          ame**2*kp4*p1p2*p1p3*p2p3**2*2D0 - 
     &          kp4*m12*p1p2*p1p3*p2p3**2*2D0 + 
     &          ame**2*p1p2*p1p3*p1p4*p2p3**2*2D0 + 
     &          ame**4*ampi**2*p1p2*p1p3*p2p4*2D0 - 
     &          ame**2*p1p2*p1p3**3*p2p4*2D0 + 
     &          ame**6*ampi**2*p2p3*p2p4*2D0 + 
     &          ame**4*ampi**2*p1p2*p2p3*p2p4*2D0 - 
     &          ame**4*p1p3**2*p2p3*p2p4*2D0 - 
     &          ame**8*ampi**2*p3p4*2D0 - 
     &          ame**6*ampi**2*p1p2*p3p4*2D0 + 
     &          ame**4*ampi**2*p1p2**2*p3p4*2D0 + 
     &          ame**2*ampi**2*p1p2**3*p3p4*2D0 + 
     &          ame**6*p1p3**2*p3p4*2D0 + 
     &          ame**2*p1p2**2*p1p3**2*p3p4*2D0 - 
     &          ame**6*p1p3*p2p3*p3p4*2D0 - 
     &          ame**2*p1p2**2*p1p3*p2p3*p3p4*2D0 + 
     &          ame**2*m12*p1p2*p1p3*p1p4*p2p3*3D0 - 
     &          ame**2*ampi**2*p1p2**2*p1p3**2*4D0 - 
     &          ame**4*p1p3**4*4D0 - 
     &          ame**6*ampi**2*p1p3*p1p4*4D0 + 
     &          ame**4*p1p3**3*p1p4*4D0 - 
     &          ame**6*ampi**2*p1p3*p2p3*4D0 + 
     &          ame**4*ampi**2*p1p2*p1p3*p2p3*4D0 + 
     &          ame**2*kp4*m12*p1p2*p1p3*p2p3*4D0 - 
     &          ame**2*m12*p1p2*p1p3**2*p2p3*4D0 - 
     &          kp4*p1p2**2*p1p3**2*p2p3*4D0 + 
     &          ame**4*p1p3**3*p2p3*4D0 - 
     &          ame**4*p1p3**2*p2p3**2*4D0 + 
     &          ame**4*p1p3*p2p3**3*4D0 - 
     &          ame**2*ampi**2*p1p2**2*p1p3*p2p4*4D0 - 
     &          ame**4*p1p3*p2p3**2*p2p4*4D0 + 
     &          ame**4*p1p2*p1p3**2*p3p4*4D0 - 
     &          ame**4*p1p2*p1p3*p2p3*p3p4*4D0 - 
     &          ame**6*ampi**2*kp4*p1p3*6D0 + 
     &          ame**2*ampi**2*kp4*p1p2**2*p1p3*6D0 + 
     &          ame**4*kp4*p1p3**3*6D0 - 
     &          ame**2*p1p2*p1p3**2*p1p4*p2p3*6D0 + 
     &          ame**4*kp4*p1p3*p2p3**2*6D0 + 
     &          ame**4*p1p3*p1p4*p2p3**2*6D0 + 
     &          ame**6*ampi**2*p1p3*p2p4*6D0 - 
     &          ame**4*p1p3**3*p2p4*6D0 + 
     &          ame**2*p1p2*p1p3**2*p2p3*p2p4*6D0 + 
     &          ame**2*ampi**2*p1p2**2*p1p3*p2p3*8D0 + 
     &          ame**2*p1p2*p1p3**3*p2p3*8D0 - 
     &          ame**2*p1p2*p1p3**2*p2p3**2*8D0 - 
     &          kp3*1D0*(-(ampi**2*p1p2**3*1D0*
     &                (m12 - p1p3*2D0)) + 
     &             p1p2*(p1p3*
     &                 (-(ame**2*p2p3*p2p4*2D0) + 
     &                   p2p3**2*(-(m12*2D0) + ame**2*4D0) - 
     &                   ame**2*p1p4*1D0*(m12 - p2p3*4D0) + 
     &                   ame**2*
     &                    (m12*p2p4*2D0 + 
     &                     ampi**2*
     &                     (-(m12*1D0) - m22*2D0 + ame**2*4D0)
     &                     )) + 
     &                p1p3**2*2D0*
     &                 (ame**2*p1p4 - ame**2*p2p4*2D0 - 
     &                   p2p3*1D0*(m22 + ame**2*6D0)) + 
     &                ame**2*
     &                 (p2p3*
     &                    (-(m12*p2p4*1D0) + 
     &                     ampi**2*
     &                     (m22 + m12*2D0 - ame**2*6D0)) + 
     &                   p1p4*
     &                    (m12*p2p3*2D0 + 
     &                     ampi**2*
     &                     (m12 + m22*2D0 - ame**2*6D0)) + 
     &                   ampi**2*
     &                    (ame**2*m12 + 
     &                     p2p4*
     &                     (-(m22*1D0) - m12*2D0 + ame**2*6D0)
     &                     ))) + 
     &             ame**2*
     &              (ame**2*ampi**2*m22*p1p4 + 
     &                ame**2*m12*p1p4*p2p3 + 
     &                m12*p2p3**2*p2p4 - 
     &                ame**2*ampi**2*m12*p2p4*1D0 + 
     &                ame**4*ampi**2*m12*2D0 - 
     &                ame**4*ampi**2*p1p4*2D0 + 
     &                ame**2*ampi**2*m12*p1p4*2D0 + 
     &                ame**2*ampi**2*m22*p2p3*2D0 - 
     &                m12*p1p4*p2p3**2*2D0 - 
     &                ame**2*ampi**2*m22*p2p4*2D0 - 
     &                ame**2*m12*p2p3*p2p4*2D0 + 
     &                p1p3**3*2D0*(m22 + ame**2*2D0) + 
     &                p1p3**2*
     &                 (m12*p2p3*2D0 + 
     &                   p2p4*2D0*(m22 - ame**2*2D0) + 
     &                   p1p4*(-(m22*1D0) + ame**2*2D0)) - 
     &                ame**4*ampi**2*p2p3*4D0 + 
     &                ame**4*ampi**2*p2p4*4D0 - 
     &                p1p3*1D0*
     &                 (ame**2*
     &                    (-(m12*p2p4*1D0) + 
     &                     ampi**2*2D0*
     &                     (m12 + m22 + ame**2*2D0)) + 
     &                   p2p3*
     &                    (ame**2*m12*2D0 + 
     &                     p2p4*
     &                     (-(m22*1D0) + ame**2*2D0 - m12*2D0)
     &                     ) + 
     &                   p1p4*
     &                    (ame**2*m12*2D0 + 
     &                     p2p3*(m12 + m22*2D0 - ame**2*4D0))
     &                    - ame**2*p2p3**2*8D0)) + 
     &             p1p2**2*
     &              (-(p1p3**2*p2p3*4D0) + 
     &                ampi**2*
     &                 (ame**2*(p2p4 - m12*1D0)*2D0 + 
     &                   p2p3*(m12 - ame**2*2D0) - 
     &                   ame**2*p1p4*4D0) + 
     &                p1p3*
     &                 (m12*p2p3*2D0 + 
     &                   ampi**2*(m22 + ame**2*1.D1)))) + 
     &          kp2*(p1p3**3*2D0*
     &              (m22*p2p3 + ame**2*p1p4*2D0 + 
     &                ame**2*(m12 - p2p4*2D0)) + 
     &             ampi**2*p1p2**2*
     &              (m12*p2p4 - m12*p1p4*1D0 - m12*p2p3*1D0 - 
     &                ame**2*ampi**2*2D0 - p1p3**2*2D0 + 
     &                ame**2*p3p4*2D0 + 
     &                p1p3*
     &                 (m12 + p1p4*2D0 + p2p3*2D0 - p2p4*2D0))
     &               - ame**2*p1p3**4*4D0 + 
     &             p1p3*(ame**2*
     &                 (ame**2*m12*(p3p4 - ampi**2*2D0) + 
     &                   p2p3*
     &                    (-(ampi**2*(m12 + m22)*2D0) + 
     &                     p3p4*(m22 - ame**2*2D0 + m12*3D0))
     &                    + p2p3**3*4D0 - p2p3**2*p2p4*4D0 + 
     &                   ampi**2*p2p4*
     &                    (-(m12*1D0) + ame**2*4D0)) + 
     &                p1p4*
     &                 (ame**2*ampi**2*(m12 - ame**2*4D0) + 
     &                   p2p3**2*(-(m12*2D0) + ame**2*4D0)))
     &              + ame**2*
     &              (-(m12*p2p3**3*2D0) + 
     &                ampi**2*p1p4*
     &                 (ame**2*m12 + p2p3*(m22 - ame**2*2D0))
     &                 + p2p3**2*
     &                 (m12*p3p4 + m12*p2p4*2D0 + 
     &                   ampi**2*
     &                    (-(m12*1D0) - m22*1D0 + ame**2*2D0))
     &                  + p2p3*
     &                 (ampi**2*p2p4*(m22 - ame**2*2D0) + 
     &                   ame**2*m12*
     &                    (-(p3p4*3D0) + ampi**2*4D0)) - 
     &                ame**2*ampi**2*1D0*
     &                 (m12*p2p4 + 
     &                   (ampi**2 - p3p4*1D0)*
     &                    (-(m12*1D0) - m22*3D0 + ame**2*6D0))
     &                ) + 
     &             p1p2*(ampi**2*p2p3**2*
     &                 (-(m12*1D0) + ame**2*2D0) + 
     &                p2p3*
     &                 (ame**2*m12*(ampi**2 - p3p4*1D0) + 
     &                   ampi**2*p1p4*(m12 - ame**2*2D0) + 
     &                   ampi**2*p2p4*(m12 - ame**2*2D0)) + 
     &                p1p3**3*p2p3*4D0 + 
     &                p1p3*
     &                 (-(ampi**2*m22*p2p4*1D0) + 
     &                   ame**2*ampi**2*p2p4*2D0 + 
     &                   p1p4*
     &                    (m12*p2p3*2D0 + 
     &                     ampi**2*(m22 - ame**2*2D0)) - 
     &                   ame**2*ampi**2*m12*3D0 + 
     &                   ame**2*m12*p3p4*3D0 + 
     &                   m12*p2p3**2*4D0 + 
     &                   p2p3*
     &                    (-(ame**2*p3p4*2D0) + 
     &                     ampi**2*
     &                     (m22 - m12*1D0 + ame**2*2D0) - 
     &                     m12*p2p4*4D0)) - 
     &                ame**2*ampi**2*1D0*(ampi**2 - p3p4*1D0)*
     &                 (-(m22*1D0) - m12*3D0 + ame**2*8D0) + 
     &                p1p3**2*
     &                 (-(p2p3*2D0*
     &                     (m12 + p1p4*2D0 - p2p4*4D0)) - 
     &                   ame**2*p3p4*6D0 - p2p3**2*8D0 + 
     &                   ampi**2*(-(m22*1D0) + ame**2*8D0)))
     &              + p1p3**2*
     &              (-(ame**2*m12*p2p3*2D0) - 
     &                p1p4*2D0*
     &                 (ame**2*m12 + p2p3*(m22 - ame**2*2D0))
     &                 + p2p3**2*(m12*2D0 - ame**2*4D0) + 
     &                ame**2*
     &                 (m12*p2p4*2D0 + 
     &                   p3p4*(m22*3D0 - ame**2*6D0) + 
     &                   ampi**2*
     &                    (-(m12*1D0) - m22*3D0 + ame**2*1.D1)
     &                   ))) - 
     &          ame**2*kp4*p1p2*p1p3**2*p2p3*1.2D1) + 
     &       kp2*(ame**2*ampi**4*m12*p1p2**2 + 
     &          ame**2*ampi**4*m22*p1p2**2 + 
     &          ampi**2*kp4*m12*p1p2**3 + 
     &          ame**4*ampi**2*kp4*m22*p1p3 + 
     &          ame**4*ampi**2*m12*p1p3**2 + 
     &          ame**2*kp4*m12*p1p2*p1p3**2 + 
     &          ame**6*ampi**2*m12*p1p4 + 
     &          ame**4*ampi**2*kp4*m12*p2p3 + 
     &          ame**2*ampi**2*m22*p1p2*p1p4*p2p3 + 
     &          ame**2*m12*p1p2*p1p3*p1p4*p2p3 + 
     &          ame**2*m22*p1p3**2*p1p4*p2p3 + 
     &          ame**4*ampi**2*m12*p2p3**2 + 
     &          ame**4*ampi**2*m22*p2p3**2 + 
     &          ame**2*kp4*m12*p1p2*p2p3**2 + 
     &          ame**2*m12*p1p2*p1p4*p2p3**2 + 
     &          ame**2*m12*p1p3*p1p4*p2p3**2 + 
     &          ame**2*m22*p1p3*p1p4*p2p3**2 + 
     &          ame**2*m12*p1p4*p2p3**3 + 
     &          ame**2*ampi**2*m12*p1p2**2*p2p4 + 
     &          ame**4*ampi**2*m22*p1p3*p2p4 + 
     &          ame**2*ampi**2*m12*p1p2*p1p3*p2p4 + 
     &          ame**4*m12*p1p3**2*p2p4 + 
     &          ame**4*m12*p1p3*p2p3*p2p4 + 
     &          ame**2*m12*p1p2*p1p3*p2p3*p2p4 + 
     &          ame**2*m22*p1p3**2*p2p3*p2p4 + 
     &          ame**2*m12*p1p3*p2p3**2*p2p4 + 
     &          ame**6*ampi**2*m12*p3p4 + 
     &          ame**6*ampi**2*m22*p3p4 + 
     &          ame**2*m12*p1p2**2*p1p3*p3p4 + 
     &          ame**2*m22*p1p2*p1p3**2*p3p4 + 
     &          ame**6*m12*p2p3*p3p4 + 
     &          ame**2*m12*p1p2*p1p3*p2p3*p3p4 - 
     &          ame**6*ampi**4*m12*1D0 - 
     &          ame**6*ampi**4*m22*1D0 - 
     &          ame**4*ampi**2*kp4*m12*p1p2*1D0 - 
     &          ampi**2*kp4*m22*p1p2**2*p1p3*1D0 - 
     &          ame**2*kp4*m22*p1p3**3*1D0 - 
     &          ame**2*ampi**2*m12*p1p2**2*p1p4*1D0 - 
     &          ame**4*ampi**2*m12*p1p3*p1p4*1D0 - 
     &          ame**2*ampi**2*m22*p1p2*p1p3*p1p4*1D0 - 
     &          ampi**2*kp4*m12*p1p2**2*p2p3*1D0 - 
     &          ame**2*kp4*m12*p1p3**2*p2p3*1D0 - 
     &          ame**4*ampi**2*m12*p1p4*p2p3*1D0 - 
     &          ame**4*m12*p1p3*p1p4*p2p3*1D0 - 
     &          ame**2*kp4*m22*p1p3*p2p3**2*1D0 - 
     &          ame**4*m12*p1p4*p2p3**2*1D0 - 
     &          ame**2*kp4*m12*p2p3**3*1D0 - 
     &          ame**6*ampi**2*m12*p2p4*1D0 - 
     &          ame**2*m12*p1p2*p1p3**2*p2p4*1D0 - 
     &          ame**2*m22*p1p3**3*p2p4*1D0 - 
     &          ame**4*ampi**2*m22*p2p3*p2p4*1D0 - 
     &          ame**2*ampi**2*m12*p1p2*p2p3*p2p4*1D0 - 
     &          ame**2*m12*p1p3**2*p2p3*p2p4*1D0 - 
     &          ame**2*ampi**2*m12*p1p2**2*p3p4*1D0 - 
     &          ame**2*ampi**2*m22*p1p2**2*p3p4*1D0 - 
     &          ame**6*m12*p1p3*p3p4*1D0 - 
     &          ame**4*m22*p1p3**2*p3p4*1D0 - 
     &          ame**2*m12*p1p2**2*p2p3*p3p4*1D0 - 
     &          ame**4*m12*p1p3*p2p3*p3p4*1D0 - 
     &          ame**4*m22*p1p3*p2p3*p3p4*1D0 - 
     &          ame**2*m22*p1p2*p1p3*p2p3*p3p4*1D0 - 
     &          ame**4*m12*p2p3**2*p3p4*1D0 - 
     &          ame**2*m12*p1p2*p2p3**2*p3p4*1D0 + 
     &          ame**8*ampi**4*2D0 - 
     &          ame**6*ampi**2*kp4*m12*2D0 + 
     &          ame**6*ampi**4*p1p2*2D0 - 
     &          ame**4*ampi**4*m12*p1p2*2D0 - 
     &          ame**4*ampi**4*p1p2**2*2D0 + 
     &          ame**2*ampi**2*kp4*m12*p1p2**2*2D0 - 
     &          ame**2*ampi**4*p1p2**3*2D0 + 
     &          ampi**4*m12*p1p2**3*2D0 + 
     &          ame**4*ampi**2*kp4*p1p2*p1p3*2D0 + 
     &          ame**4*ampi**2*m12*p1p2*p1p3*2D0 - 
     &          ampi**2*kp4*p1p2**3*p1p3*2D0 - 
     &          ampi**2*m12*p1p2**3*p1p3*2D0 - 
     &          ame**6*ampi**2*p1p3**2*2D0 + 
     &          ame**4*kp4*m12*p1p3**2*2D0 + 
     &          ame**2*ampi**2*m12*p1p2*p1p3**2*2D0 - 
     &          ampi**2*m22*p1p2**2*p1p3**2*2D0 - 
     &          ame**2*kp4*p1p2*p1p3**3*2D0 - 
     &          ame**2*m12*p1p2*p1p3**3*2D0 - 
     &          ame**2*m22*p1p3**4*2D0 - 
     &          ame**4*ampi**2*m22*p1p3*p1p4*2D0 - 
     &          ame**2*ampi**2*m12*p1p2*p1p3*p1p4*2D0 + 
     &          ame**2*m12*p1p2*p1p3**2*p1p4*2D0 + 
     &          ame**2*m22*p1p3**3*p1p4*2D0 + 
     &          ame**4*ampi**2*m12*p1p3*p2p3*2D0 - 
     &          ame**2*ampi**2*m12*p1p2*p1p3*p2p3*2D0 - 
     &          ame**2*ampi**2*m22*p1p2*p1p3*p2p3*2D0 - 
     &          kp4*m12*p1p2**2*p1p3*p2p3*2D0 + 
     &          kp4*m22*p1p2*p1p3**2*p2p3*2D0 - 
     &          ame**2*m12*p1p3**3*p2p3*2D0 - 
     &          ame**4*ampi**2*p1p2*p1p4*p2p3*2D0 - 
     &          ame**2*ampi**2*p1p2**2*p1p4*p2p3*2D0 + 
     &          ampi**2*m12*p1p2**2*p1p4*p2p3*2D0 - 
     &          m12*p1p2**2*p1p3*p1p4*p2p3*2D0 - 
     &          ame**4*p1p3**2*p1p4*p2p3*2D0 + 
     &          ame**2*m12*p1p3**2*p1p4*p2p3*2D0 - 
     &          m22*p1p2*p1p3**2*p1p4*p2p3*2D0 - 
     &          ame**6*ampi**2*p2p3**2*2D0 + 
     &          ame**4*kp4*m12*p2p3**2*2D0 - 
     &          ame**4*ampi**2*p1p2*p2p3**2*2D0 + 
     &          ame**2*ampi**2*m12*p1p2*p2p3**2*2D0 - 
     &          ame**2*kp4*p1p2*p1p3*p2p3**2*2D0 - 
     &          ame**2*m12*p1p2*p1p3*p2p3**2*2D0 + 
     &          kp4*m12*p1p2*p1p3*p2p3**2*2D0 - 
     &          ame**2*m22*p1p3**2*p2p3**2*2D0 + 
     &          ame**4*p1p3*p1p4*p2p3**2*2D0 - 
     &          ame**2*p1p2*p1p3*p1p4*p2p3**2*2D0 - 
     &          m12*p1p2*p1p3*p1p4*p2p3**2*2D0 - 
     &          ame**2*m12*p1p3*p2p3**3*2D0 + 
     &          ame**6*ampi**2*p1p3*p2p4*2D0 + 
     &          ame**2*ampi**2*m22*p1p2*p1p3*p2p4*2D0 + 
     &          ampi**2*m12*p1p2**2*p1p3*p2p4*2D0 - 
     &          m12*p1p2**2*p1p3**2*p2p4*2D0 - 
     &          ame**4*p1p3**3*p2p4*2D0 - 
     &          m22*p1p2*p1p3**3*p2p4*2D0 + 
     &          ame**6*ampi**2*p2p3*p2p4*2D0 + 
     &          ame**4*ampi**2*p1p2*p2p3*p2p4*2D0 - 
     &          ame**4*p1p3**2*p2p3*p2p4*2D0 - 
     &          m12*p1p2*p1p3**2*p2p3*p2p4*2D0 - 
     &          ame**8*ampi**2*p3p4*2D0 - 
     &          ame**6*ampi**2*p1p2*p3p4*2D0 + 
     &          ame**4*ampi**2*m12*p1p2*p3p4*2D0 + 
     &          ame**4*ampi**2*p1p2**2*p3p4*2D0 + 
     &          ame**2*ampi**2*p1p2**3*p3p4*2D0 - 
     &          ampi**2*m12*p1p2**3*p3p4*2D0 - 
     &          ame**4*m12*p1p2*p1p3*p3p4*2D0 + 
     &          m12*p1p2**3*p1p3*p3p4*2D0 + 
     &          ame**6*p1p3**2*p3p4*2D0 + 
     &          m22*p1p2**2*p1p3**2*p3p4*2D0 + 
     &          ame**6*p1p3*p2p3*p3p4*2D0 + 
     &          ame**2*p1p2**2*p1p3*p2p3*p3p4*2D0 + 
     &          m12*p1p2**2*p1p3*p2p3*p3p4*2D0 + 
     &          ame**4*ampi**2*m22*p1p3**2*3D0 + 
     &          ampi**2*p1p2**3*p1p3**2*4D0 + 
     &          ame**2*p1p2*p1p3**4*4D0 - 
     &          ame**2*p1p2*p1p3**3*p1p4*4D0 - 
     &          ame**6*ampi**2*p1p3*p2p3*4D0 + 
     &          ame**4*ampi**2*p1p2*p1p3*p2p3*4D0 - 
     &          ame**2*kp4*m12*p1p2*p1p3*p2p3*4D0 + 
     &          kp4*p1p2**2*p1p3**2*p2p3*4D0 + 
     &          m12*p1p2**2*p1p3**2*p2p3*4D0 + 
     &          ame**4*p1p3**3*p2p3*4D0 + 
     &          m22*p1p2*p1p3**3*p2p3*4D0 + 
     &          p1p2**2*p1p3**2*p1p4*p2p3*4D0 - 
     &          ame**2*p1p2*p1p3**2*p2p3**2*4D0 + 
     &          m12*p1p2*p1p3**2*p2p3**2*4D0 + 
     &          ame**4*p1p3*p2p3**3*4D0 + 
     &          p1p2**2*p1p3**3*p2p4*4D0 - 
     &          ame**4*p1p3*p2p3**2*p2p4*4D0 - 
     &          p1p2**3*p1p3**2*p3p4*4D0 + 
     &          ame**4*p1p2*p1p3*p2p3*p3p4*4D0 + 
     &          ame**6*ampi**2*kp4*p1p3*6D0 - 
     &          ame**2*ampi**2*kp4*p1p2**2*p1p3*6D0 - 
     &          ame**4*ampi**2*p1p2*p1p3**2*6D0 - 
     &          ame**4*kp4*p1p3**3*6D0 + 
     &          ame**4*ampi**2*p1p2*p1p3*p1p4*6D0 + 
     &          ame**2*ampi**2*p1p2**2*p1p3*p1p4*6D0 - 
     &          ampi**2*m12*p1p2**2*p1p3*p2p3*6D0 - 
     &          ame**2*p1p2*p1p3**2*p1p4*p2p3*6D0 - 
     &          ame**4*kp4*p1p3*p2p3**2*6D0 - 
     &          ame**4*ampi**2*p1p2*p1p3*p2p4*6D0 + 
     &          ame**2*p1p2*p1p3**3*p2p4*6D0 + 
     &          ame**2*p1p2*p1p3**2*p2p3*p2p4*6D0 - 
     &          ame**2*p1p2**2*p1p3**2*p3p4*6D0 + 
     &          kp3*(-(ampi**2*p1p2**3*1D0*
     &                (m12 - p1p3*2D0)) + 
     &             ame**2*
     &              (p1p3**3*2D0*(m22 + ame**2*2D0) + 
     &                m12*(ame**2*ampi**2 - p2p3**2*1D0)*
     &                 (p2p4 + ame**2*2D0 - p2p3*2D0) + 
     &                ame**2*p1p4*
     &                 (m12*p2p3 + ampi**2*(m22 - ame**2*2D0))
     &                  + p1p3**2*
     &                 (m12*(p2p3 - ame**2*1D0)*2D0 + 
     &                   p1p4*(-(m22*1D0) + ame**2*2D0)) - 
     &                p1p3*1D0*
     &                 (m12*p1p4*p2p3 + 
     &                   p2p3*p2p4*(m22 - ame**2*2D0) - 
     &                   p2p3**2*2D0*(m22 + ame**2*2D0) + 
     &                   ame**2*
     &                    (m12*p2p4 + 
     &                     ampi**2*2D0*(m22 + ame**2*2D0))))
     &              + p1p2*
     &              (p1p3**2*2D0*
     &                 (ame**2*p1p4 + 
     &                   p2p4*(m22 - ame**2*2D0) - 
     &                   p2p3*2D0*(m22 + ame**2*2D0)) + 
     &                ame**2*
     &                 (ampi**2*p1p4*(m12 - ame**2*2D0) + 
     &                   p2p3*
     &                    (-(m12*p2p4*1D0) + 
     &                     ampi**2*(m22 - ame**2*2D0)) + 
     &                   ampi**2*
     &                    (ame**2*m12 + 
     &                     p2p4*(-(m22*1D0) + ame**2*2D0))) - 
     &                p1p3*1D0*
     &                 (ame**2*ampi**2*m12 + 
     &                   ame**2*m12*p1p4 - 
     &                   p2p3*2D0*
     &                    ((ame**2 + m12)*p2p4 + 
     &                     ame**2*m12*2D0) + m12*p2p3**2*4D0))
     &               + p1p2**2*
     &              (ampi**2*
     &                 (-((ame**2*m12 + 
     &                     p2p4*(m12 - ame**2*1D0))*2D0) + 
     &                   p2p3*(-(ame**2*2D0) + m12*3D0)) - 
     &                p1p3**2*p2p4*4D0 + 
     &                p1p3*
     &                 (m12*p2p4*2D0 + 
     &                   ampi**2*(m22 + ame**2*6D0)))) + 
     &          ame**2*ampi**2*p1p2**2*p1p3*p2p3*8D0 - 
     &          p1p2**2*p1p3**3*p2p3*8D0 - 
     &          ame**2*ampi**2*p1p2**2*p1p3*p2p4*8D0 + 
     &          kp2*(ame**2*
     &              (m12*(ame**2*ampi**2 - p2p3**2*1D0)*
     &                 (ampi**2 - p3p4*1D0) + 
     &                p1p3**3*2D0*(m12 + p1p4*2D0) + 
     &                p1p4*
     &                 (ame**2*ampi**2*m12 - 
     &                   m12*p2p3**2*2D0 + 
     &                   ampi**2*p2p3*
     &                    (-(m22*1D0) + ame**2*2D0)) - 
     &                p1p3**4*4D0 + 
     &                p1p3**2*
     &                 (-(m12*p1p4*2D0) - p2p3**2*4D0 + 
     &                   ampi**2*(-(m12*1D0) + ame**2*4D0)) + 
     &                p1p3*
     &                 (m12*p2p3**2*2D0 + 
     &                   p2p3*p3p4*(m22 - ame**2*2D0) + 
     &                   ame**2*m12*(p3p4 - ampi**2*2D0) + 
     &                   p1p4*
     &                    (p2p3**2*4D0 + 
     &                     ampi**2*(m12 - ame**2*4D0)))) + 
     &             p1p2**2*
     &              (ampi**2*
     &                 (-(m12*p1p4*1D0) + 
     &                   (ame**2 - m12*1D0)*
     &                    (ampi**2 - p3p4*1D0)*2D0) + 
     &                p1p3*
     &                 (ampi**2*p1p4*2D0 + 
     &                   m12*(-(p3p4*2D0) + ampi**2*3D0)) + 
     &                p1p3**2*(p3p4*4D0 - ampi**2*6D0)) + 
     &             p1p2*(ame**2*ampi**2*(ampi**2 - p3p4*1D0)*
     &                 (-(m22*1D0) + ame**2*2D0) + 
     &                p2p3*
     &                 (ame**2*m12*(p3p4 - ampi**2*1D0) + 
     &                   ampi**2*p1p4*
     &                    (-(m12*1D0) + ame**2*2D0)) + 
     &                p1p3**2*
     &                 (-(1D0*(-(m22*1D0) + ame**2*2D0)*
     &                     (ampi**2 - p3p4*2D0)) - 
     &                   p2p3*(m12 + p1p4*2D0)*4D0) + 
     &                p1p3*
     &                 (p2p3*
     &                    (-((ame**2 + m12)*p3p4*2D0) + 
     &                     ampi**2*m12*3D0) + 
     &                   p1p4*
     &                    (ampi**2*(m22 - ame**2*2D0) + 
     &                     m12*p2p3*4D0)) + p1p3**3*p2p3*8D0))
     &            + ame**2*kp4*p1p2*p1p3**2*p2p3*1.2D1)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*
     &     (ampi**2*p1p2**2 + 
     &       ame**2*(p1p3**2 + p2p3**2 - 
     &          ame**2*ampi**2*1D0) - p1p2*p1p3*p2p3*2D0))
