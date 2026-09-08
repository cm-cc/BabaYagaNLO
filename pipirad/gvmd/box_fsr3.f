        elmat2 = (ampi**2*D3aex12x35x14epem1m2p*2D0*
     &     (kp4*(p1p3*p2p4 - kp2*p1p3*1D0 + 
     &          ame**2*(ampi**2 - p3p4*1D0) + 
     &          p1p2*(ampi**2 - p3p4*1D0) + 
     &          p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)) + 
     &       kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &          (ame**2 + p1p2)*p3p4*1D0 + 
     &          p2p4*(p1p3 - kp1*1D0) + 
     &          p1p4*(p2p3 - kp2*1D0 - p2p4*2D0) - 
     &          kp4*2D0*(p1p2 + ame**2*3D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)) - 
     &  (C0ax12px35m2m1p*1D0*
     &     (kp4*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)) - 
     &  (D2aex12x35x14epem1m2p*
     &     (kp3*(p2p4*((ame**2 + p1p2)*p3p4 - 
     &             ampi**2*(ame**2 + p1p2)*1D0 + 
     &             kp4*(ame**2 - p1p2*1D0) + 
     &             p2p4*(kp1 - p1p3*1D0)) + 
     &          p1p4**2*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          p1p4*((ame**2 + p1p2)*p3p4 - 
     &             ampi**2*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(kp1 + kp2 - p1p3*1D0 - p2p3*1D0) + 
     &             p2p4**2*2D0 + kp4*2D0*(p1p2 + ame**2*2D0)))
     &         + kp4*(-(p1p4**2*p2p3*1D0) + 
     &          p2p3*(ampi**2*(ame**2 + p1p2) + 
     &             (kp1 + p1p3)*p2p4 - 
     &             (ame**2 + p1p2)*(kp4 + p3p4)*1D0) + 
     &          p1p4*(kp2*(p1p3 + p2p3) + 
     &             ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             p2p3*(kp1 + p2p3 + p1p3*2D0) - 
     &             p2p4*1D0*(p1p3 + p2p3*2D0))))*4D0)/
     &   (kp3*kp4*(ame**2 + p1p2)) - 
     &  (D0aex12x35x14epem1m2p*p1p4*
     &     (kp4*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0)))*4D0)/
     &   (kp3*kp4*(ame**2 + p1p2)) - 
     &  (D1aeex35px12x14m1em2p*2D0*
     &     (ame**2*ampi**2*kp4*p1p3 + 
     &       ame**2*ampi**2*kp3*p2p3 + 
     &       ame**2*ampi**2*kp4*p2p3 + 
     &       ame**2*ampi**2*kp4*p2p4 - 
     &       ame**2*ampi**2*kp3*p1p3*1D0 - 
     &       ame**2*kp4**2*p1p3*1D0 - 
     &       ame**2*ampi**2*kp3*p1p4*1D0 - 
     &       ame**2*ampi**2*kp4*p1p4*1D0 - 
     &       ame**2*kp4**2*p2p3*1D0 - 
     &       ame**2*ampi**2*kp3*p2p4*1D0 - 
     &       ampi**2*kp3*p1p2*p1p4*2D0 - 
     &       ampi**2*kp4*p1p2*p1p4*2D0 - 
     &       kp3*p1p4**2*p2p3*2D0 - kp4*p1p4**2*p2p3*2D0 - 
     &       kp3*p1p3*p1p4*p2p4*2D0 - 
     &       kp4*p1p3*p1p4*p2p4*2D0 + 
     &       p3p4*(kp3*(ame**2 + p1p2)*p1p4 + 
     &          kp4*(p1p2*p1p4 - ame**2*p2p3*1D0))*2D0 - 
     &       kp2*1D0*(kp4*(ame**2*p3p4 - p1p3*p1p4*2D0) + 
     &          kp3*(ame**2*ampi**2 - p1p4**2*2D0)) + 
     &       kp1*(kp4*(ame**2*p3p4 + p1p4*p2p3*2D0) + 
     &          kp3*(ame**2*ampi**2 + p1p4*p2p4*2D0)) - 
     &       ame**2*kp3*kp4*p2p4*3D0 + 
     &       kp3*kp4*p1p2*p1p4*4D0 + kp4*p1p3*p1p4*p2p3*4D0 + 
     &       kp3*p1p4**2*p2p4*4D0 + ame**2*kp3*kp4*p1p4*5D0))/
     &   (kp3*kp4*(ame**2 + p1p2))
                  
