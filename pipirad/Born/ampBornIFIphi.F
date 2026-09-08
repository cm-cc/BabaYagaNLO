        elmat = (2.5D-1*(Cf2*kp1**3*(p2p3 - p2p4*1D0)*
     &       (ampi**2 - p3p4*1D0) - 
     &      kp1**2*1D0*(ame**2*Cf1*kp4 + ame**2*Cf1*p1p4 + 
     &         Cf1*p1p2*p2p3 + Cf2*kp4*p1p4*p2p3 + 
     &         Cf2*kp4*p1p3*p2p4 + ame**2*Cf2*kp4*p3p4 + 
     &         Cf2*kp4*p1p2*p3p4 - 
     &         ame**2*ampi**2*Cf2*kp4*1D0 - 
     &         ampi**2*Cf2*kp4*p1p2*1D0 - 
     &         ame**2*Cf1*p1p3*1D0 - Cf2*kp4*p1p3*p2p3*1D0 - 
     &         Cf2*kp4*p2p3**2*1D0 - Cf1*p1p2*p2p4*1D0 - 
     &         Cf2*kp4*p1p4*p2p4*1D0 - Cf2*kp4*p2p4**2*1D0 + 
     &         ame**2*Cf1*p2p3*2D0 + 
     &         ame**2*ampi**2*Cf2*p2p3*2D0 + 
     &         ampi**2*Cf2*p1p2*p2p3*2D0 + 
     &         ame**2*Cf3*p2p3**2*2D0 + 
     &         Cf3*p1p2*p2p3**2*2D0 - Cf2*p1p3*p2p3**2*2D0 + 
     &         Cf2*p1p4*p2p3**2*2D0 - ame**2*Cf1*p2p4*2D0 - 
     &         ame**2*ampi**2*Cf2*p2p4*2D0 - 
     &         ampi**2*Cf2*p1p2*p2p4*2D0 + 
     &         Cf2*kp4*p2p3*p2p4*2D0 + 
     &         ame**2*Cf3*p2p4**2*2D0 + 
     &         Cf3*p1p2*p2p4**2*2D0 - Cf2*p1p3*p2p4**2*2D0 + 
     &         Cf2*p1p4*p2p4**2*2D0 - 
     &         ame**2*Cf2*p2p3*p3p4*2D0 - 
     &         Cf2*p1p2*p2p3*p3p4*2D0 + 
     &         ame**2*Cf2*p2p4*p3p4*2D0 + 
     &         Cf2*p1p2*p2p4*p3p4*2D0 + 
     &         kp3*(ame**2*ampi**2*Cf2 + Cf2*p2p3**2 + 
     &            Cf2*p1p4*p2p4 + Cf2*p2p4**2 - 
     &            ame**2*Cf1*1D0 - Cf2*p1p4*p2p3*1D0 - 
     &            ame**2*Cf2*p3p4*1D0 + 
     &            Cf2*p1p3*(p2p3 - p2p4*1D0) + 
     &            Cf2*p1p2*(ampi**2 - p3p4*1D0) - 
     &            Cf2*p2p3*p2p4*2D0) - 
     &         ame**2*Cf3*p2p3*p2p4*4D0 - 
     &         Cf3*p1p2*p2p3*p2p4*4D0 + 
     &         Cf2*p1p3*p2p3*p2p4*4D0 - 
     &         Cf2*p1p4*p2p3*p2p4*4D0 + 
     &         kp2*(ampi**2*Cf2*p1p4 + Cf2*p1p3*p3p4 - 
     &            ampi**2*Cf2*p1p3*1D0 - Cf2*p1p4*p3p4*1D0 + 
     &            Cf1*p2p3*2D0 - Cf1*p2p4*2D0 + 
     &            Cf2*kp4*(p3p4 - ampi**2*1D0)*2D0 + 
     &            Cf2*kp3*(ampi**2 - p3p4*1D0)*2D0 + 
     &            ame**2*ampi**2*Cf3*4D0 + 
     &            ampi**2*Cf3*p1p2*4D0 - 
     &            ame**2*Cf3*p3p4*4D0 - Cf3*p1p2*p3p4*4D0)) + 
     &      kp1*((kp3 - kp4*1D0)*
     &          (ame**2*Cf2*kp3**2 + ame**2*Cf2*kp4**2 - 
     &            kp4*1D0*
     &             (ame**2*
     &                (ame**2*Cf3 + Cf2*p1p3 - Cf2*p1p4*1D0)*
     &                2D0 + 
     &               p1p2*
     &                (Cf2*p1p3 + Cf2*p2p3 - Cf2*p1p4*1D0 - 
     &                  Cf2*p2p4*1D0 + ame**2*Cf3*2D0)) + 
     &            kp3*(-(ame**2*Cf2*kp4*2D0) + 
     &               ame**2*
     &                (ame**2*Cf3 + Cf2*p1p3 - Cf2*p1p4*1D0)*
     &                2D0 + 
     &               p1p2*
     &                (Cf2*p1p3 + Cf2*p2p3 - Cf2*p1p4*1D0 - 
     &                  Cf2*p2p4*1D0 + ame**2*Cf3*2D0)) + 
     &            (ame**2 + p1p2)*
     &             (ame**2*Cf1 + ame**2*ampi**2*Cf2*2D0 - 
     &               ame**2*Cf3*p1p4*2D0 + 
     &               Cf2*p1p4*p2p3*2D0 - Cf2*p1p4*p2p4*2D0 - 
     &               ame**2*Cf2*p3p4*2D0 + 
     &               p1p3*
     &                (ame**2*Cf3 + Cf2*p2p4 - Cf2*p2p3*1D0)*
     &                2D0 + 
     &               p1p2*
     &                (Cf1 + ampi**2*Cf2*2D0 + Cf3*p2p3*2D0 - 
     &                  Cf3*p2p4*2D0 - Cf2*p3p4*2D0))) + 
     &         kp2**2*(ampi**2*Cf2*p2p4 + Cf2*p2p3*p3p4 - 
     &            ampi**2*Cf2*p2p3*1D0 - Cf2*p2p4*p3p4*1D0 + 
     &            Cf1*p1p3*2D0 - Cf1*p1p4*2D0 + 
     &            Cf2*kp4*(p3p4 - ampi**2*1D0)*2D0 + 
     &            Cf2*kp3*(ampi**2 - p3p4*1D0)*2D0 + 
     &            ame**2*ampi**2*Cf3*4D0 + 
     &            ampi**2*Cf3*p1p2*4D0 - 
     &            ame**2*Cf3*p3p4*4D0 - Cf3*p1p2*p3p4*4D0) + 
     &         kp2*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &          (ame**2*Cf1 - Cf1*p1p2*1D0 + 
     &            ame**2*ampi**2*Cf2*2D0 + Cf2*kp3**2*2D0 + 
     &            Cf2*kp4**2*2D0 + ampi**2*Cf2*p1p2*2D0 - 
     &            Cf2*p1p3*p2p3*2D0 + Cf2*p1p4*p2p3*2D0 + 
     &            Cf2*p1p3*p2p4*2D0 - Cf2*p1p4*p2p4*2D0 - 
     &            ame**2*Cf2*p3p4*2D0 - Cf2*p1p2*p3p4*2D0 - 
     &            kp4*1D0*
     &             (Cf2*p1p4 + Cf2*p2p4 - Cf2*p1p3*1D0 - 
     &               Cf2*p2p3*1D0 + ame**2*Cf3*4D0 + 
     &               Cf3*p1p2*4D0) + 
     &            kp3*(Cf2*p1p4 + Cf2*p2p4 - Cf2*p1p3*1D0 - 
     &               Cf2*p2p3*1D0 + ame**2*Cf3*4D0 - 
     &               Cf2*kp4*4D0 + Cf3*p1p2*4D0))) + 
     &      kp2*(ame**2*Cf2*kp4**3 - ame**2*Cf2*kp3**3*1D0 - 
     &         kp4**2*1D0*
     &          (ame**2*(ame**2*Cf3 + Cf2*p2p3 - 
     &               Cf2*p2p4*1D0)*2D0 + 
     &            p1p2*(Cf2*p1p3 + Cf2*p2p3 - Cf2*p1p4*1D0 - 
     &               Cf2*p2p4*1D0 + ame**2*Cf3*2D0)) + 
     &         kp4*(kp2*(ame**2*Cf1 + Cf2*p1p4*p2p3 + 
     &               ame**2*Cf2*p3p4 - 
     &               ame**2*ampi**2*Cf2*1D0 - 
     &               Cf2*p1p3**2*1D0 - Cf2*p1p4**2*1D0 - 
     &               Cf2*p1p4*p2p4*1D0 + 
     &               Cf2*p1p2*(p3p4 - ampi**2*1D0) + 
     &               Cf2*p1p3*(p2p4 - p2p3*1D0 + p1p4*2D0)) + 
     &            (ame**2 + p1p2)*
     &             (ame**2*Cf1 + ame**2*ampi**2*Cf2*2D0 - 
     &               ame**2*Cf3*p2p4*2D0 + 
     &               Cf2*p1p3*p2p4*2D0 - Cf2*p1p4*p2p4*2D0 - 
     &               ame**2*Cf2*p3p4*2D0 + 
     &               p2p3*
     &                (ame**2*Cf3 + Cf2*p1p4 - Cf2*p1p3*1D0)*
     &                2D0 + 
     &               p1p2*
     &                (Cf1 + ampi**2*Cf2*2D0 + Cf3*p1p3*2D0 - 
     &                  Cf3*p1p4*2D0 - Cf2*p3p4*2D0))) + 
     &         kp3**2*(-(ame**2*
     &               (ame**2*Cf3 + Cf2*p2p3 - Cf2*p2p4*1D0)*
     &               2D0) + 
     &            p1p2*(Cf2*p1p4 + Cf2*p2p4 - Cf2*p1p3*1D0 - 
     &               Cf2*p2p3*1D0 - ame**2*Cf3*2D0) + 
     &            ame**2*Cf2*kp4*3D0) - 
     &         kp3*1D0*(-(kp4*2D0*
     &               (ame**2*
     &                  (ame**2*Cf3 + Cf2*p2p3 - 
     &                    Cf2*p2p4*1D0)*2D0 + 
     &                 p1p2*
     &                  (Cf2*p1p3 + Cf2*p2p3 - Cf2*p1p4*1D0 - 
     &                    Cf2*p2p4*1D0 + ame**2*Cf3*2D0))) + 
     &            kp2*(ame**2*Cf1 + Cf2*p1p4*p2p3 + 
     &               ame**2*Cf2*p3p4 - 
     &               ame**2*ampi**2*Cf2*1D0 - 
     &               Cf2*p1p3**2*1D0 - Cf2*p1p4**2*1D0 - 
     &               Cf2*p1p4*p2p4*1D0 + 
     &               Cf2*p1p2*(p3p4 - ampi**2*1D0) + 
     &               Cf2*p1p3*(p2p4 - p2p3*1D0 + p1p4*2D0)) + 
     &            (ame**2 + p1p2)*
     &             (ame**2*Cf1 + ame**2*ampi**2*Cf2*2D0 - 
     &               ame**2*Cf3*p2p4*2D0 + 
     &               Cf2*p1p3*p2p4*2D0 - Cf2*p1p4*p2p4*2D0 - 
     &               ame**2*Cf2*p3p4*2D0 + 
     &               p2p3*
     &                (ame**2*Cf3 + Cf2*p1p4 - Cf2*p1p3*1D0)*
     &                2D0 + 
     &               p1p2*
     &                (Cf1 + ampi**2*Cf2*2D0 + Cf3*p1p3*2D0 - 
     &                  Cf3*p1p4*2D0 - Cf2*p3p4*2D0)) + 
     &            ame**2*Cf2*kp4**2*3D0) + 
     &         kp2*(ame**2*Cf1*(p2p4 - p2p3*1D0) + 
     &            p1p3**2*
     &             (ame**2*Cf3 + Cf3*p1p2 + Cf2*p2p4 - 
     &               Cf2*p2p3*1D0)*2D0 + 
     &            p1p4**2*
     &             (ame**2*Cf3 + Cf3*p1p2 + Cf2*p2p4 - 
     &               Cf2*p2p3*1D0)*2D0 - 
     &            p1p4*1D0*
     &             (Cf2*kp2*(p3p4 - ampi**2*1D0) + 
     &               ame**2*
     &                (Cf1 + ampi**2*Cf2 - Cf2*p3p4*1D0)*2D0
     &                + p1p2*
     &                (Cf1 + ampi**2*Cf2*2D0 - Cf2*p3p4*2D0))
     &             + p1p3*
     &             (Cf2*kp2*(p3p4 - ampi**2*1D0) + 
     &               ame**2*
     &                (Cf1 + ampi**2*Cf2 - Cf2*p3p4*1D0)*2D0
     &                - p1p4*
     &                (ame**2*Cf3 + Cf2*p2p4 - Cf2*p2p3*1D0)*
     &                4D0 + 
     &               p1p2*
     &                (Cf1 + ampi**2*Cf2*2D0 - Cf2*p3p4*2D0 - 
     &                  Cf3*p1p4*4D0))))))/
     &  (kp1*kp2*(ame**2 + p1p2)*(ampi**2 + p3p4))
                  
