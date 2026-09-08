        elmat = (2.5D-1*((kp3 - kp4*1D0)**2*2D0*
     &       (ame**2*Cf2*f2*kp3**2 + ame**2*Cf2*f2*kp4**2 - 
     &         kp4*(ame**2 + p1p2)*1D0*
     &          (Cf2*f2*p1p3 + Cf2*f2*p2p3 - 
     &            Cf2*f2*p1p4*1D0 - Cf2*f2*p2p4*1D0 + 
     &            ame**2*Cf3*f2*2D0 + ame**2*Cf2*f3*2D0) + 
     &         kp3*(-(ame**2*Cf2*f2*kp4*2D0) + 
     &            (ame**2 + p1p2)*
     &             (Cf2*f2*p1p3 + Cf2*f2*p2p3 - 
     &               Cf2*f2*p1p4*1D0 - Cf2*f2*p2p4*1D0 + 
     &               ame**2*Cf3*f2*2D0 + ame**2*Cf2*f3*2D0))
     &          + (ame**2 + p1p2)*
     &          (ame**2*Cf2*f1 + ame**2*Cf1*f2 + 
     &            ame**2*Cf3*f2*p2p3 + ame**2*Cf2*f3*p2p3 - 
     &            ame**2*Cf3*f2*p1p4*1D0 - 
     &            ame**2*Cf2*f3*p1p4*1D0 - 
     &            ame**2*Cf3*f2*p2p4*1D0 - 
     &            ame**2*Cf2*f3*p2p4*1D0 + 
     &            ame**2*ampi**2*Cf2*f2*2D0 + 
     &            Cf2*f2*p1p4*p2p3*2D0 - 
     &            Cf2*f2*p1p4*p2p4*2D0 - 
     &            ame**2*Cf2*f2*p3p4*2D0 + 
     &            p1p3*(ame**2*(Cf3*f2 + Cf2*f3) - 
     &               Cf2*f2*p2p3*2D0 + Cf2*f2*p2p4*2D0) + 
     &            ame**4*Cf3*f3*4D0 + 
     &            p1p2*(Cf2*f1 + Cf1*f2 + 
     &               (Cf3*f2 + Cf2*f3)*p1p3 + Cf3*f2*p2p3 + 
     &               Cf2*f3*p2p3 - 
     &               (Cf3*f2 + Cf2*f3)*p1p4*1D0 - 
     &               Cf3*f2*p2p4*1D0 - Cf2*f3*p2p4*1D0 + 
     &               ampi**2*Cf2*f2*2D0 - Cf2*f2*p3p4*2D0 + 
     &               ame**2*Cf3*f3*4D0))) + 
     &      kp1**2*((Cf2*f1 + Cf1*f2)*p2p3**2 + 
     &         Cf2*f1*kp4*p2p4 + Cf1*f2*kp4*p2p4 + 
     &         Cf2*f1*p1p3*p2p4 + Cf1*f2*p1p3*p2p4 + 
     &         Cf2*f1*p2p4**2 + Cf1*f2*p2p4**2 - 
     &         Cf2*f1*kp3*p2p4*1D0 - Cf1*f2*kp3*p2p4*1D0 - 
     &         Cf2*f1*p1p4*p2p4*1D0 - Cf1*f2*p1p4*p2p4*1D0 + 
     &         ame**2*ampi**2*Cf2*f1*2D0 + 
     &         ame**2*ampi**2*Cf1*f2*2D0 - 
     &         ame**2*Cf3*f1*p2p4*2D0 - 
     &         ame**2*Cf1*f3*p2p4*2D0 - 
     &         ame**2*Cf2*f1*p3p4*2D0 - 
     &         ame**2*Cf1*f2*p3p4*2D0 + ame**2*Cf1*f1*4D0 + 
     &         ame**2*ampi**4*Cf2*f2*4D0 - 
     &         ame**2*ampi**2*Cf3*f2*p2p4*4D0 - 
     &         ame**2*ampi**2*Cf2*f3*p2p4*4D0 - 
     &         ampi**2*Cf2*f2*kp3*p2p4*4D0 + 
     &         ampi**2*Cf2*f2*kp4*p2p4*4D0 + 
     &         ampi**2*Cf2*f2*p1p3*p2p4*4D0 - 
     &         ampi**2*Cf2*f2*p1p4*p2p4*4D0 + 
     &         ame**2*Cf3*f2*p2p4*p3p4*4D0 + 
     &         ame**2*Cf2*f3*p2p4*p3p4*4D0 + 
     &         Cf2*f2*kp3*p2p4*p3p4*4D0 - 
     &         Cf2*f2*kp4*p2p4*p3p4*4D0 - 
     &         Cf2*f2*p1p3*p2p4*p3p4*4D0 + 
     &         Cf2*f2*p1p4*p2p4*p3p4*4D0 + 
     &         ame**2*Cf2*f2*p3p4**2*4D0 + 
     &         p1p2*2D0*(Cf1*f1 + ampi**2*Cf2*f1 + 
     &            ampi**2*Cf1*f2 - Cf2*f1*p3p4*1D0 - 
     &            Cf1*f2*p3p4*1D0 + ampi**4*Cf2*f2*2D0 + 
     &            Cf2*f2*p3p4**2*2D0 + 
     &            p2p3*(Cf3*f1 + Cf1*f3 + 
     &               ampi**2*Cf3*f2*2D0 + 
     &               ampi**2*Cf2*f3*2D0 - 
     &               (Cf3*f2 + Cf2*f3)*p3p4*2D0) - 
     &            p2p4*1D0*
     &             (Cf3*f1 + Cf1*f3 + ampi**2*Cf3*f2*2D0 + 
     &               ampi**2*Cf2*f3*2D0 - 
     &               (Cf3*f2 + Cf2*f3)*p3p4*2D0) - 
     &            ampi**2*Cf2*f2*p3p4*4D0) + 
     &         p2p3*(Cf2*f1*p1p4 + Cf1*f2*p1p4 - 
     &            Cf2*f1*p1p3*1D0 - Cf1*f2*p1p3*1D0 + 
     &            ame**2*Cf3*f1*2D0 + ame**2*Cf1*f3*2D0 - 
     &            Cf2*f1*p2p4*2D0 - Cf1*f2*p2p4*2D0 + 
     &            ame**2*ampi**2*Cf3*f2*4D0 + 
     &            ame**2*ampi**2*Cf2*f3*4D0 - 
     &            ampi**2*Cf2*f2*p1p3*4D0 + 
     &            ampi**2*Cf2*f2*p1p4*4D0 - 
     &            ame**2*Cf3*f2*p3p4*4D0 - 
     &            ame**2*Cf2*f3*p3p4*4D0 + 
     &            Cf2*f2*p1p3*p3p4*4D0 - 
     &            Cf2*f2*p1p4*p3p4*4D0 + 
     &            kp3*(Cf2*f1 + Cf1*f2 + ampi**2*Cf2*f2*4D0 - 
     &               Cf2*f2*p3p4*4D0) - 
     &            kp4*1D0*
     &             (Cf2*f1 + Cf1*f2 + ampi**2*Cf2*f2*4D0 - 
     &               Cf2*f2*p3p4*4D0)) - 
     &         ame**2*ampi**2*Cf2*f2*p3p4*8D0) + 
     &      kp2**2*((Cf2*f1 + Cf1*f2)*p1p3**2 + 
     &         Cf2*f1*kp4*p1p4 + Cf1*f2*kp4*p1p4 + 
     &         Cf2*f1*p1p4**2 + Cf1*f2*p1p4**2 + 
     &         Cf2*f1*p1p4*p2p3 + Cf1*f2*p1p4*p2p3 - 
     &         Cf2*f1*kp3*p1p4*1D0 - Cf1*f2*kp3*p1p4*1D0 - 
     &         Cf2*f1*p1p4*p2p4*1D0 - Cf1*f2*p1p4*p2p4*1D0 + 
     &         ame**2*ampi**2*Cf2*f1*2D0 + 
     &         ame**2*ampi**2*Cf1*f2*2D0 - 
     &         ame**2*Cf3*f1*p1p4*2D0 - 
     &         ame**2*Cf1*f3*p1p4*2D0 - 
     &         ame**2*Cf2*f1*p3p4*2D0 - 
     &         ame**2*Cf1*f2*p3p4*2D0 + ame**2*Cf1*f1*4D0 + 
     &         ame**2*ampi**4*Cf2*f2*4D0 - 
     &         ame**2*ampi**2*Cf3*f2*p1p4*4D0 - 
     &         ame**2*ampi**2*Cf2*f3*p1p4*4D0 - 
     &         ampi**2*Cf2*f2*kp3*p1p4*4D0 + 
     &         ampi**2*Cf2*f2*kp4*p1p4*4D0 + 
     &         ampi**2*Cf2*f2*p1p4*p2p3*4D0 - 
     &         ampi**2*Cf2*f2*p1p4*p2p4*4D0 + 
     &         ame**2*Cf3*f2*p1p4*p3p4*4D0 + 
     &         ame**2*Cf2*f3*p1p4*p3p4*4D0 + 
     &         Cf2*f2*kp3*p1p4*p3p4*4D0 - 
     &         Cf2*f2*kp4*p1p4*p3p4*4D0 - 
     &         Cf2*f2*p1p4*p2p3*p3p4*4D0 + 
     &         Cf2*f2*p1p4*p2p4*p3p4*4D0 + 
     &         ame**2*Cf2*f2*p3p4**2*4D0 + 
     &         p1p2*2D0*(Cf1*f1 + ampi**2*Cf2*f1 + 
     &            ampi**2*Cf1*f2 - Cf2*f1*p3p4*1D0 - 
     &            Cf1*f2*p3p4*1D0 + ampi**4*Cf2*f2*2D0 + 
     &            Cf2*f2*p3p4**2*2D0 + 
     &            p1p3*(Cf3*f1 + Cf1*f3 + 
     &               ampi**2*Cf3*f2*2D0 + 
     &               ampi**2*Cf2*f3*2D0 - 
     &               (Cf3*f2 + Cf2*f3)*p3p4*2D0) - 
     &            p1p4*1D0*
     &             (Cf3*f1 + Cf1*f3 + ampi**2*Cf3*f2*2D0 + 
     &               ampi**2*Cf2*f3*2D0 - 
     &               (Cf3*f2 + Cf2*f3)*p3p4*2D0) - 
     &            ampi**2*Cf2*f2*p3p4*4D0) + 
     &         p1p3*(Cf2*f1*p2p4 + Cf1*f2*p2p4 - 
     &            Cf2*f1*p2p3*1D0 - Cf1*f2*p2p3*1D0 + 
     &            ame**2*Cf3*f1*2D0 + ame**2*Cf1*f3*2D0 - 
     &            Cf2*f1*p1p4*2D0 - Cf1*f2*p1p4*2D0 + 
     &            ame**2*ampi**2*Cf3*f2*4D0 + 
     &            ame**2*ampi**2*Cf2*f3*4D0 - 
     &            ampi**2*Cf2*f2*p2p3*4D0 + 
     &            ampi**2*Cf2*f2*p2p4*4D0 - 
     &            ame**2*Cf3*f2*p3p4*4D0 - 
     &            ame**2*Cf2*f3*p3p4*4D0 + 
     &            Cf2*f2*p2p3*p3p4*4D0 - 
     &            Cf2*f2*p2p4*p3p4*4D0 + 
     &            kp3*(Cf2*f1 + Cf1*f2 + ampi**2*Cf2*f2*4D0 - 
     &               Cf2*f2*p3p4*4D0) - 
     &            kp4*1D0*
     &             (Cf2*f1 + Cf1*f2 + ampi**2*Cf2*f2*4D0 - 
     &               Cf2*f2*p3p4*4D0)) - 
     &         ame**2*ampi**2*Cf2*f2*p3p4*8D0) + 
     &      kp2*(kp3 - kp4*1D0)*
     &       (ame**2*Cf2*f1*p2p4 + ame**2*Cf1*f2*p2p4 + 
     &         Cf2*f1*p1p2*p2p4 + Cf1*f2*p1p2*p2p4 - 
     &         ame**2*Cf2*f1*p2p3*1D0 - 
     &         ame**2*Cf1*f2*p2p3*1D0 - 
     &         Cf2*f1*p1p2*p2p3*1D0 - Cf1*f2*p1p2*p2p3*1D0 + 
     &         Cf3*f1*p1p2**2*2D0 + Cf1*f3*p1p2**2*2D0 - 
     &         ame**2*Cf3*f2*p1p3**2*2D0 - 
     &         ame**2*Cf2*f3*p1p3**2*2D0 - 
     &         Cf3*f2*p1p2*p1p3**2*2D0 - 
     &         Cf2*f3*p1p2*p1p3**2*2D0 - 
     &         ame**2*Cf3*f2*p1p4**2*2D0 - 
     &         ame**2*Cf2*f3*p1p4**2*2D0 - 
     &         Cf3*f2*p1p2*p1p4**2*2D0 - 
     &         Cf2*f3*p1p2*p1p4**2*2D0 + 
     &         Cf2*f2*kp3**2*(p1p3 - p1p4*1D0)*2D0 + 
     &         Cf2*f2*kp4**2*(p1p3 - p1p4*1D0)*2D0 - 
     &         ame**2*Cf2*f1*p1p3*3D0 - 
     &         ame**2*Cf1*f2*p1p3*3D0 - 
     &         Cf2*f1*p1p2*p1p3*3D0 - Cf1*f2*p1p2*p1p3*3D0 + 
     &         ame**2*Cf2*f1*p1p4*3D0 + 
     &         ame**2*Cf1*f2*p1p4*3D0 + 
     &         Cf2*f1*p1p2*p1p4*3D0 + Cf1*f2*p1p2*p1p4*3D0 - 
     &         ame**2*ampi**2*Cf2*f2*p1p3*4D0 - 
     &         ampi**2*Cf2*f2*p1p2*p1p3*4D0 + 
     &         ame**2*ampi**2*Cf2*f2*p1p4*4D0 + 
     &         ampi**2*Cf2*f2*p1p2*p1p4*4D0 + 
     &         ame**2*Cf3*f2*p1p3*p1p4*4D0 + 
     &         ame**2*Cf2*f3*p1p3*p1p4*4D0 + 
     &         Cf3*f2*p1p2*p1p3*p1p4*4D0 + 
     &         Cf2*f3*p1p2*p1p3*p1p4*4D0 - 
     &         ame**2*ampi**2*Cf2*f2*p2p3*4D0 - 
     &         ampi**2*Cf2*f2*p1p2*p2p3*4D0 + 
     &         Cf2*f2*p1p3**2*p2p3*4D0 + 
     &         Cf2*f2*p1p4**2*p2p3*4D0 + 
     &         Cf2*f2*p1p3*p2p3**2*4D0 - 
     &         Cf2*f2*p1p4*p2p3**2*4D0 + 
     &         ame**2*ampi**2*Cf2*f2*p2p4*4D0 + 
     &         ampi**2*Cf2*f2*p1p2*p2p4*4D0 - 
     &         Cf2*f2*p1p3**2*p2p4*4D0 - 
     &         Cf2*f2*p1p4**2*p2p4*4D0 + 
     &         Cf2*f2*p1p3*p2p4**2*4D0 - 
     &         Cf2*f2*p1p4*p2p4**2*4D0 + 
     &         ame**2*Cf2*f2*p1p3*p3p4*4D0 + 
     &         Cf2*f2*p1p2*p1p3*p3p4*4D0 - 
     &         ame**2*Cf2*f2*p1p4*p3p4*4D0 - 
     &         Cf2*f2*p1p2*p1p4*p3p4*4D0 + 
     &         ame**2*Cf2*f2*p2p3*p3p4*4D0 + 
     &         Cf2*f2*p1p2*p2p3*p3p4*4D0 - 
     &         ame**2*Cf2*f2*p2p4*p3p4*4D0 - 
     &         Cf2*f2*p1p2*p2p4*p3p4*4D0 + 
     &         ame**4*Cf3*f1*6D0 + ame**4*Cf1*f3*6D0 - 
     &         ame**2*Cf3*f2*p1p3*p2p3*6D0 - 
     &         ame**2*Cf2*f3*p1p3*p2p3*6D0 - 
     &         Cf3*f2*p1p2*p1p3*p2p3*6D0 - 
     &         Cf2*f3*p1p2*p1p3*p2p3*6D0 + 
     &         ame**2*Cf3*f2*p1p4*p2p3*6D0 + 
     &         ame**2*Cf2*f3*p1p4*p2p3*6D0 + 
     &         Cf3*f2*p1p2*p1p4*p2p3*6D0 + 
     &         Cf2*f3*p1p2*p1p4*p2p3*6D0 + 
     &         ame**2*Cf3*f2*p1p3*p2p4*6D0 + 
     &         ame**2*Cf2*f3*p1p3*p2p4*6D0 + 
     &         Cf3*f2*p1p2*p1p3*p2p4*6D0 + 
     &         Cf2*f3*p1p2*p1p3*p2p4*6D0 - 
     &         ame**2*Cf3*f2*p1p4*p2p4*6D0 - 
     &         ame**2*Cf2*f3*p1p4*p2p4*6D0 - 
     &         Cf3*f2*p1p2*p1p4*p2p4*6D0 - 
     &         Cf2*f3*p1p2*p1p4*p2p4*6D0 - 
     &         kp4*1D0*(-(Cf2*f2*p1p3**2*2D0) - 
     &            Cf2*f2*p1p4**2*2D0 + ame**2*Cf2*f1*3D0 + 
     &            ame**2*Cf1*f2*3D0 - 
     &            ame**2*Cf3*f2*p1p4*4D0 - 
     &            ame**2*Cf2*f3*p1p4*4D0 + 
     &            p1p2*(Cf2*f1 + Cf1*f2 + 
     &               (Cf3*f2 + Cf2*f3)*p1p3*4D0 - 
     &               (Cf3*f2 + Cf2*f3)*p1p4*4D0) + 
     &            Cf2*f2*p1p4*p2p3*6D0 - 
     &            Cf2*f2*p1p4*p2p4*6D0 + 
     &            p1p3*(ame**2*Cf3*f2*4D0 + 
     &               ame**2*Cf2*f3*4D0 + Cf2*f2*p1p4*4D0 - 
     &               Cf2*f2*p2p3*6D0 + Cf2*f2*p2p4*6D0)) + 
     &         kp3*(-(Cf2*f2*p1p3**2*2D0) - 
     &            Cf2*f2*p1p4**2*2D0 + ame**2*Cf2*f1*3D0 + 
     &            ame**2*Cf1*f2*3D0 - 
     &            ame**2*Cf3*f2*p1p4*4D0 - 
     &            ame**2*Cf2*f3*p1p4*4D0 + 
     &            Cf2*f2*kp4*p1p4*4D0 + 
     &            p1p2*(Cf2*f1 + Cf1*f2 + 
     &               (Cf3*f2 + Cf2*f3)*p1p3*4D0 - 
     &               (Cf3*f2 + Cf2*f3)*p1p4*4D0) + 
     &            Cf2*f2*p1p4*p2p3*6D0 - 
     &            Cf2*f2*p1p4*p2p4*6D0 + 
     &            p1p3*(ame**2*Cf3*f2*4D0 + 
     &               ame**2*Cf2*f3*4D0 - Cf2*f2*kp4*4D0 + 
     &               Cf2*f2*p1p4*4D0 - Cf2*f2*p2p3*6D0 + 
     &               Cf2*f2*p2p4*6D0)) + 
     &         ame**2*Cf3*f1*p1p2*8D0 + 
     &         ame**2*Cf1*f3*p1p2*8D0 + 
     &         ame**4*Cf3*f3*p1p3*8D0 + 
     &         Cf3*f3*p1p2**2*p1p3*8D0 - 
     &         ame**4*Cf3*f3*p1p4*8D0 - 
     &         Cf3*f3*p1p2**2*p1p4*8D0 - 
     &         Cf2*f2*p1p3*p1p4*p2p3*8D0 + 
     &         Cf2*f2*p1p3*p1p4*p2p4*8D0 - 
     &         Cf2*f2*p1p3*p2p3*p2p4*8D0 + 
     &         Cf2*f2*p1p4*p2p3*p2p4*8D0 + 
     &         ame**2*Cf3*f3*p1p2*p1p3*1.6D1 - 
     &         ame**2*Cf3*f3*p1p2*p1p4*1.6D1) + 
     &      kp1*((kp3 - kp4*1D0)*
     &          (ame**2*Cf2*f1*p1p4 + ame**2*Cf1*f2*p1p4 + 
     &            Cf2*f1*p1p2*p1p4 + Cf1*f2*p1p2*p1p4 - 
     &            ame**2*Cf2*f1*p1p3*1D0 - 
     &            ame**2*Cf1*f2*p1p3*1D0 - 
     &            Cf2*f1*p1p2*p1p3*1D0 - 
     &            Cf1*f2*p1p2*p1p3*1D0 + Cf3*f1*p1p2**2*2D0 + 
     &            Cf1*f3*p1p2**2*2D0 - 
     &            ame**2*Cf3*f2*p2p3**2*2D0 - 
     &            ame**2*Cf2*f3*p2p3**2*2D0 - 
     &            Cf3*f2*p1p2*p2p3**2*2D0 - 
     &            Cf2*f3*p1p2*p2p3**2*2D0 - 
     &            ame**2*Cf3*f2*p2p4**2*2D0 - 
     &            ame**2*Cf2*f3*p2p4**2*2D0 - 
     &            Cf3*f2*p1p2*p2p4**2*2D0 - 
     &            Cf2*f3*p1p2*p2p4**2*2D0 + 
     &            Cf2*f2*kp3**2*(p2p3 - p2p4*1D0)*2D0 + 
     &            Cf2*f2*kp4**2*(p2p3 - p2p4*1D0)*2D0 - 
     &            ame**2*Cf2*f1*p2p3*3D0 - 
     &            ame**2*Cf1*f2*p2p3*3D0 - 
     &            Cf2*f1*p1p2*p2p3*3D0 - 
     &            Cf1*f2*p1p2*p2p3*3D0 + 
     &            ame**2*Cf2*f1*p2p4*3D0 + 
     &            ame**2*Cf1*f2*p2p4*3D0 + 
     &            Cf2*f1*p1p2*p2p4*3D0 + 
     &            Cf1*f2*p1p2*p2p4*3D0 - 
     &            ame**2*ampi**2*Cf2*f2*p1p3*4D0 - 
     &            ampi**2*Cf2*f2*p1p2*p1p3*4D0 + 
     &            ame**2*ampi**2*Cf2*f2*p1p4*4D0 + 
     &            ampi**2*Cf2*f2*p1p2*p1p4*4D0 - 
     &            ame**2*ampi**2*Cf2*f2*p2p3*4D0 - 
     &            ampi**2*Cf2*f2*p1p2*p2p3*4D0 + 
     &            Cf2*f2*p1p3**2*p2p3*4D0 + 
     &            Cf2*f2*p1p4**2*p2p3*4D0 + 
     &            Cf2*f2*p1p3*p2p3**2*4D0 - 
     &            Cf2*f2*p1p4*p2p3**2*4D0 + 
     &            ame**2*ampi**2*Cf2*f2*p2p4*4D0 + 
     &            ampi**2*Cf2*f2*p1p2*p2p4*4D0 - 
     &            Cf2*f2*p1p3**2*p2p4*4D0 - 
     &            Cf2*f2*p1p4**2*p2p4*4D0 + 
     &            ame**2*Cf3*f2*p2p3*p2p4*4D0 + 
     &            ame**2*Cf2*f3*p2p3*p2p4*4D0 + 
     &            Cf3*f2*p1p2*p2p3*p2p4*4D0 + 
     &            Cf2*f3*p1p2*p2p3*p2p4*4D0 + 
     &            Cf2*f2*p1p3*p2p4**2*4D0 - 
     &            Cf2*f2*p1p4*p2p4**2*4D0 + 
     &            ame**2*Cf2*f2*p1p3*p3p4*4D0 + 
     &            Cf2*f2*p1p2*p1p3*p3p4*4D0 - 
     &            ame**2*Cf2*f2*p1p4*p3p4*4D0 - 
     &            Cf2*f2*p1p2*p1p4*p3p4*4D0 + 
     &            ame**2*Cf2*f2*p2p3*p3p4*4D0 + 
     &            Cf2*f2*p1p2*p2p3*p3p4*4D0 - 
     &            ame**2*Cf2*f2*p2p4*p3p4*4D0 - 
     &            Cf2*f2*p1p2*p2p4*p3p4*4D0 + 
     &            ame**4*Cf3*f1*6D0 + ame**4*Cf1*f3*6D0 - 
     &            ame**2*Cf3*f2*p1p3*p2p3*6D0 - 
     &            ame**2*Cf2*f3*p1p3*p2p3*6D0 - 
     &            Cf3*f2*p1p2*p1p3*p2p3*6D0 - 
     &            Cf2*f3*p1p2*p1p3*p2p3*6D0 + 
     &            ame**2*Cf3*f2*p1p4*p2p3*6D0 + 
     &            ame**2*Cf2*f3*p1p4*p2p3*6D0 + 
     &            Cf3*f2*p1p2*p1p4*p2p3*6D0 + 
     &            Cf2*f3*p1p2*p1p4*p2p3*6D0 + 
     &            ame**2*Cf3*f2*p1p3*p2p4*6D0 + 
     &            ame**2*Cf2*f3*p1p3*p2p4*6D0 + 
     &            Cf3*f2*p1p2*p1p3*p2p4*6D0 + 
     &            Cf2*f3*p1p2*p1p3*p2p4*6D0 - 
     &            ame**2*Cf3*f2*p1p4*p2p4*6D0 - 
     &            ame**2*Cf2*f3*p1p4*p2p4*6D0 - 
     &            Cf3*f2*p1p2*p1p4*p2p4*6D0 - 
     &            Cf2*f3*p1p2*p1p4*p2p4*6D0 - 
     &            kp4*1D0*
     &             (-(Cf2*f2*p2p3**2*2D0) - 
     &               Cf2*f2*p2p4**2*2D0 + ame**2*Cf2*f1*3D0 + 
     &               ame**2*Cf1*f2*3D0 - 
     &               ame**2*Cf3*f2*p2p4*4D0 - 
     &               ame**2*Cf2*f3*p2p4*4D0 + 
     &               p1p2*
     &                (Cf2*f1 + Cf1*f2 + 
     &                  (Cf3*f2 + Cf2*f3)*p2p3*4D0 - 
     &                  (Cf3*f2 + Cf2*f3)*p2p4*4D0) + 
     &               Cf2*f2*p1p3*p2p4*6D0 - 
     &               Cf2*f2*p1p4*p2p4*6D0 + 
     &               p2p3*
     &                (ame**2*(Cf3*f2 + Cf2*f3)*4D0 + 
     &                  Cf2*f2*p2p4*4D0 - Cf2*f2*p1p3*6D0 + 
     &                  Cf2*f2*p1p4*6D0)) + 
     &            kp3*(-(Cf2*f2*p2p3**2*2D0) - 
     &               Cf2*f2*p2p4**2*2D0 + ame**2*Cf2*f1*3D0 + 
     &               ame**2*Cf1*f2*3D0 - 
     &               ame**2*Cf3*f2*p2p4*4D0 - 
     &               ame**2*Cf2*f3*p2p4*4D0 + 
     &               Cf2*f2*kp4*p2p4*4D0 + 
     &               p1p2*
     &                (Cf2*f1 + Cf1*f2 + 
     &                  (Cf3*f2 + Cf2*f3)*p2p3*4D0 - 
     &                  (Cf3*f2 + Cf2*f3)*p2p4*4D0) + 
     &               Cf2*f2*p1p3*p2p4*6D0 - 
     &               Cf2*f2*p1p4*p2p4*6D0 + 
     &               p2p3*
     &                (ame**2*Cf3*f2*4D0 + 
     &                  ame**2*Cf2*f3*4D0 - Cf2*f2*kp4*4D0 + 
     &                  Cf2*f2*p2p4*4D0 - Cf2*f2*p1p3*6D0 + 
     &                  Cf2*f2*p1p4*6D0)) + 
     &            ame**2*Cf3*f1*p1p2*8D0 + 
     &            ame**2*Cf1*f3*p1p2*8D0 + 
     &            ame**4*Cf3*f3*p2p3*8D0 + 
     &            Cf3*f3*p1p2**2*p2p3*8D0 - 
     &            Cf2*f2*p1p3*p1p4*p2p3*8D0 - 
     &            ame**4*Cf3*f3*p2p4*8D0 - 
     &            Cf3*f3*p1p2**2*p2p4*8D0 + 
     &            Cf2*f2*p1p3*p1p4*p2p4*8D0 - 
     &            Cf2*f2*p1p3*p2p3*p2p4*8D0 + 
     &            Cf2*f2*p1p4*p2p3*p2p4*8D0 + 
     &            ame**2*Cf3*f3*p1p2*p2p3*1.6D1 - 
     &            ame**2*Cf3*f3*p1p2*p2p4*1.6D1) + 
     &         kp2*(Cf2*f1*p1p3**2 + Cf1*f2*p1p3**2 + 
     &            Cf2*f1*p1p4**2 + Cf1*f2*p1p4**2 + 
     &            Cf2*f1*p2p3**2 + Cf1*f2*p2p3**2 + 
     &            Cf2*f1*p2p4**2 + Cf1*f2*p2p4**2 - 
     &            ame**2*Cf3*f1*p1p3*2D0 - 
     &            ame**2*Cf1*f3*p1p3*2D0 - 
     &            Cf3*f1*p1p2*p1p3*2D0 - 
     &            Cf1*f3*p1p2*p1p3*2D0 + 
     &            ame**2*Cf3*f1*p1p4*2D0 + 
     &            ame**2*Cf1*f3*p1p4*2D0 + 
     &            Cf3*f1*p1p2*p1p4*2D0 + 
     &            Cf1*f3*p1p2*p1p4*2D0 - 
     &            Cf2*f1*p1p3*p1p4*2D0 - 
     &            Cf1*f2*p1p3*p1p4*2D0 - 
     &            ame**2*Cf3*f1*p2p3*2D0 - 
     &            ame**2*Cf1*f3*p2p3*2D0 - 
     &            Cf3*f1*p1p2*p2p3*2D0 - 
     &            Cf1*f3*p1p2*p2p3*2D0 - 
     &            Cf2*f1*p1p3*p2p3*2D0 - 
     &            Cf1*f2*p1p3*p2p3*2D0 + 
     &            Cf2*f1*p1p4*p2p3*2D0 + 
     &            Cf1*f2*p1p4*p2p3*2D0 + 
     &            ame**2*Cf3*f1*p2p4*2D0 + 
     &            ame**2*Cf1*f3*p2p4*2D0 + 
     &            Cf3*f1*p1p2*p2p4*2D0 + 
     &            Cf1*f3*p1p2*p2p4*2D0 + 
     &            Cf2*f1*p1p3*p2p4*2D0 + 
     &            Cf1*f2*p1p3*p2p4*2D0 - 
     &            Cf2*f1*p1p4*p2p4*2D0 - 
     &            Cf1*f2*p1p4*p2p4*2D0 - 
     &            Cf2*f1*p2p3*p2p4*2D0 - 
     &            Cf1*f2*p2p3*p2p4*2D0 + ame**2*Cf1*f1*4D0 + 
     &            ame**2*ampi**2*Cf2*f1*4D0 + 
     &            ame**2*ampi**2*Cf1*f2*4D0 + 
     &            ampi**2*Cf2*f1*p1p2*4D0 + 
     &            ampi**2*Cf1*f2*p1p2*4D0 + 
     &            ame**2*ampi**2*Cf3*f2*p1p3*4D0 + 
     &            ame**2*ampi**2*Cf2*f3*p1p3*4D0 + 
     &            ampi**2*Cf3*f2*p1p2*p1p3*4D0 + 
     &            ampi**2*Cf2*f3*p1p2*p1p3*4D0 - 
     &            ame**2*ampi**2*Cf3*f2*p1p4*4D0 - 
     &            ame**2*ampi**2*Cf2*f3*p1p4*4D0 - 
     &            ampi**2*Cf3*f2*p1p2*p1p4*4D0 - 
     &            ampi**2*Cf2*f3*p1p2*p1p4*4D0 + 
     &            ame**2*ampi**2*Cf3*f2*p2p3*4D0 + 
     &            ame**2*ampi**2*Cf2*f3*p2p3*4D0 + 
     &            ampi**2*Cf3*f2*p1p2*p2p3*4D0 + 
     &            ampi**2*Cf2*f3*p1p2*p2p3*4D0 - 
     &            ame**2*ampi**2*Cf3*f2*p2p4*4D0 - 
     &            ame**2*ampi**2*Cf2*f3*p2p4*4D0 - 
     &            ampi**2*Cf3*f2*p1p2*p2p4*4D0 - 
     &            ampi**2*Cf2*f3*p1p2*p2p4*4D0 - 
     &            ame**2*Cf2*f1*p3p4*4D0 - 
     &            ame**2*Cf1*f2*p3p4*4D0 - 
     &            Cf2*f1*p1p2*p3p4*4D0 - 
     &            Cf1*f2*p1p2*p3p4*4D0 - 
     &            ame**2*Cf3*f2*p1p3*p3p4*4D0 - 
     &            ame**2*Cf2*f3*p1p3*p3p4*4D0 - 
     &            Cf3*f2*p1p2*p1p3*p3p4*4D0 - 
     &            Cf2*f3*p1p2*p1p3*p3p4*4D0 + 
     &            ame**2*Cf3*f2*p1p4*p3p4*4D0 + 
     &            ame**2*Cf2*f3*p1p4*p3p4*4D0 + 
     &            Cf3*f2*p1p2*p1p4*p3p4*4D0 + 
     &            Cf2*f3*p1p2*p1p4*p3p4*4D0 - 
     &            ame**2*Cf3*f2*p2p3*p3p4*4D0 - 
     &            ame**2*Cf2*f3*p2p3*p3p4*4D0 - 
     &            Cf3*f2*p1p2*p2p3*p3p4*4D0 - 
     &            Cf2*f3*p1p2*p2p3*p3p4*4D0 + 
     &            ame**2*Cf3*f2*p2p4*p3p4*4D0 + 
     &            ame**2*Cf2*f3*p2p4*p3p4*4D0 + 
     &            Cf3*f2*p1p2*p2p4*p3p4*4D0 + 
     &            Cf2*f3*p1p2*p2p4*p3p4*4D0 + 
     &            Cf2*f2*kp3**2*(p3p4 - ampi**2*1D0)*4D0 + 
     &            Cf2*f2*kp4**2*(p3p4 - ampi**2*1D0)*4D0 + 
     &            ame**2*ampi**4*Cf2*f2*8D0 + 
     &            ampi**4*Cf2*f2*p1p2*8D0 - 
     &            ampi**2*Cf2*f2*p1p3*p2p3*8D0 + 
     &            ampi**2*Cf2*f2*p1p4*p2p3*8D0 + 
     &            ampi**2*Cf2*f2*p1p3*p2p4*8D0 - 
     &            ampi**2*Cf2*f2*p1p4*p2p4*8D0 + 
     &            Cf2*f2*p1p3*p2p3*p3p4*8D0 - 
     &            Cf2*f2*p1p4*p2p3*p3p4*8D0 - 
     &            Cf2*f2*p1p3*p2p4*p3p4*8D0 + 
     &            Cf2*f2*p1p4*p2p4*p3p4*8D0 + 
     &            ame**2*Cf2*f2*p3p4**2*8D0 + 
     &            Cf2*f2*p1p2*p3p4**2*8D0 + 
     &            kp3*(Cf2*f1*p1p4 + Cf1*f2*p1p4 + 
     &               Cf2*f1*p2p4 + Cf1*f2*p2p4 - 
     &               Cf2*f1*p1p3*1D0 - Cf1*f2*p1p3*1D0 - 
     &               Cf2*f1*p2p3*1D0 - Cf1*f2*p2p3*1D0 + 
     &               ampi**2*Cf2*f2*p1p3*4D0 - 
     &               ampi**2*Cf2*f2*p1p4*4D0 + 
     &               ampi**2*Cf2*f2*p2p3*4D0 - 
     &               ampi**2*Cf2*f2*p2p4*4D0 - 
     &               Cf2*f2*p1p3*p3p4*4D0 + 
     &               Cf2*f2*p1p4*p3p4*4D0 - 
     &               Cf2*f2*p2p3*p3p4*4D0 + 
     &               Cf2*f2*p2p4*p3p4*4D0 - 
     &               ame**2*ampi**2*Cf3*f2*8D0 - 
     &               ame**2*ampi**2*Cf2*f3*8D0 + 
     &               ame**2*Cf3*f2*p3p4*8D0 + 
     &               ame**2*Cf2*f3*p3p4*8D0 + 
     &               Cf2*f2*kp4*(ampi**2 - p3p4*1D0)*8D0 - 
     &               (Cf3*f2 + Cf2*f3)*p1p2*
     &                (ampi**2 - p3p4*1D0)*8D0) + 
     &            kp4*(Cf2*f1*p2p3 + Cf1*f2*p2p3 - 
     &               Cf2*f1*p1p4*1D0 - Cf1*f2*p1p4*1D0 - 
     &               Cf2*f1*p2p4*1D0 - Cf1*f2*p2p4*1D0 + 
     &               ampi**2*Cf2*f2*p1p4*4D0 - 
     &               ampi**2*Cf2*f2*p2p3*4D0 + 
     &               ampi**2*Cf2*f2*p2p4*4D0 - 
     &               Cf2*f2*p1p4*p3p4*4D0 + 
     &               Cf2*f2*p2p3*p3p4*4D0 - 
     &               Cf2*f2*p2p4*p3p4*4D0 + 
     &               p1p3*
     &                (Cf2*f1 + Cf1*f2 - ampi**2*Cf2*f2*4D0 + 
     &                  Cf2*f2*p3p4*4D0) + 
     &               ame**2*ampi**2*Cf3*f2*8D0 + 
     &               ame**2*ampi**2*Cf2*f3*8D0 - 
     &               ame**2*Cf3*f2*p3p4*8D0 - 
     &               ame**2*Cf2*f3*p3p4*8D0 + 
     &               (Cf3*f2 + Cf2*f3)*p1p2*
     &                (ampi**2 - p3p4*1D0)*8D0) - 
     &            ame**4*ampi**2*Cf3*f3*1.6D1 - 
     &            ampi**2*Cf3*f3*p1p2**2*1.6D1 - 
     &            ame**2*ampi**2*Cf2*f2*p3p4*1.6D1 + 
     &            ame**4*Cf3*f3*p3p4*1.6D1 - 
     &            ampi**2*Cf2*f2*p1p2*p3p4*1.6D1 + 
     &            Cf3*f3*p1p2**2*p3p4*1.6D1 - 
     &            ame**2*ampi**2*Cf3*f3*p1p2*3.2D1 + 
     &            ame**2*Cf3*f3*p1p2*p3p4*3.2D1))))/
     &  (ame**2 + p1p2)**2
