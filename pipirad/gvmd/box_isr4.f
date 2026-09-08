        elmat1 = (C0ax12px35m1m2p*
     &     (kp1**2*(p2p4 - p2p3*1D0) + 
     &       kp2*(ame**2*p1p3 + kp2*p1p3 + ame**2*p2p4 - 
     &          ame**2*p1p4*1D0 - kp2*p1p4*1D0 - 
     &          ame**2*p2p3*1D0 + kp4*(p1p2 + ame**2*2D0) - 
     &          kp3*1D0*(p1p2 + ame**2*2D0)) + 
     &       kp1*((ame**2 + kp2)*
     &           (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0) + 
     &          kp3*(p1p2 + ame**2*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*2D0))))/
     &   (kp1*kp2*(ampi**2 + p3p4)) + 
     &  (ampi**2*D3aex12px24ex35em1m2p*2D0*
     &     (kp1**2*(p2p4 - p2p3*1D0) + 
     &       kp2*(ame**2*p1p3 + kp2*p1p3 + ame**2*p2p4 - 
     &          ame**2*p1p4*1D0 - kp2*p1p4*1D0 - 
     &          ame**2*p2p3*1D0 + kp4*(p1p2 + ame**2*2D0) - 
     &          kp3*1D0*(p1p2 + ame**2*2D0)) + 
     &       kp1*((ame**2 + kp2)*
     &           (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0) + 
     &          kp3*(p1p2 + ame**2*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*2D0))))/
     &   (kp1*kp2*(ampi**2 + p3p4)) + 
     &  (D1aeepx35x12x24m1em2p*2D0*
     &     (kp1**2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(p2p4 - p2p3*1D0)*2D0) + 
     &       kp2*(kp2*(p2p4*(p1p3 - p1p4*1D0) + 
     &             ame**2*(ampi**2 - p3p4*1D0))*2D0 + 
     &          p1p2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &             p2p4*(kp4 - kp3*1D0)*2D0) + 
     &          ame**2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p4*(kp3 - p2p3*1D0) + p2p4**2*2D0 - 
     &             kp4*1D0*(p1p3 - p2p3*2D0) + 
     &             p2p4*(p1p3 - kp3*2D0 - p2p3*2D0))) + 
     &       kp1*(p1p2*(ame**2*(ampi**2 - p3p4*1D0) + 
     &             p2p4*(kp3 - kp4*1D0)*2D0) + 
     &          kp2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &             p2p4*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &              2D0) + 
     &          ame**2*(p1p4*(p2p3 - kp3*1D0) + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p2p4*(p1p3 + kp3*2D0 - p1p4*2D0) - 
     &             kp4*1D0*(p1p3 - p1p4*2D0 + p2p4*2D0)))))/
     &   (kp1*kp2*(ampi**2 + p3p4)) + 
     &  (D0aex12px24ex35em1m2p*p2p4*
     &     (kp1**2*(p2p4 - p2p3*1D0) + 
     &       kp2*(ame**2*p1p3 + kp2*p1p3 + ame**2*p2p4 - 
     &          ame**2*p1p4*1D0 - kp2*p1p4*1D0 - 
     &          ame**2*p2p3*1D0 + kp4*(p1p2 + ame**2*2D0) - 
     &          kp3*1D0*(p1p2 + ame**2*2D0)) + 
     &       kp1*((ame**2 + kp2)*
     &           (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0) + 
     &          kp3*(p1p2 + ame**2*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*2D0)))*4D0)/
     &   (kp1*kp2*(ampi**2 + p3p4)) + 
     &  (D1aex12px24ex35em1m2p*2D0*
     &     (kp1**2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0)*2D0 - 
     &          (p1p4 + p2p4)*(p2p3 - p2p4*1D0)*2D0) + 
     &       kp2*(p1p2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &             p2p4*(kp4 - kp3*1D0))*2D0 + 
     &          kp2*(ame**2*(ampi**2 - p3p4*1D0) + 
     &             p2p4*(p1p3 - p1p4*1D0)*2D0) + 
     &          ame**2*(ame**2*(p3p4 - ampi**2*1D0)*2D0 - 
     &             (p1p4 + p2p4)*(p2p3 - p2p4*1D0)*2D0 + 
     &             kp4*(p1p3 + p2p3 + p2p4*2D0) - 
     &             kp3*1D0*(p1p4 + p2p4*3D0))) + 
     &       kp1*(p1p2**2*(ampi**2 - p3p4*1D0)*2D0 + 
     &          kp2*(p2p4**2 + p1p3*(p1p4 + p2p4) - 
     &             p1p4**2*1D0 - (p1p4 + p2p3)*p2p4*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0))*2D0 + 
     &          p1p2*2D0*
     &           ((kp3 + p1p3)*p2p4 - 
     &             kp4*1D0*(p1p3 + p2p4 - p1p4*1D0) + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p4*(p2p3 - p2p4*2D0)) + 
     &          ame**2*(p2p4*
     &              (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*2D0 + 
     &             kp3*(p1p4 + p2p4*3D0) + 
     &             kp4*(p2p3 - p1p3*1D0 - p2p4*4D0)))))/
     &   (kp1*kp2*(ampi**2 + p3p4))
