        elmat1 = -((C12apx450pm1p*
     &       (kp2*(kp4 - p1p4*1D0)*
     &          ((ame**2 + p1p2)*(kp3 - kp4*1D0) + 
     &            kp2*(p1p4 - p1p3*1D0)) + 
     &         kp1**2*((kp4 - p2p4*1D0)*(p2p3 - p2p4*1D0) + 
     &            kp2*(ampi**2 - p3p4*1D0)) + 
     &         kp1*(kp2**2*(p3p4 - ampi**2*1D0) - 
     &            (ame**2 + p1p2)*1D0*(kp3 - kp4*1D0)*
     &             (kp4 - p2p4*1D0) + 
     &            kp2*(-(p2p4*1D0*(kp3 + p1p3 - kp4*1D0)) + 
     &               p1p4*(kp3 + p2p3 - kp4*1D0))))*2D0)/
     &     (kp1*kp2*(ampi**2 + p3p4)*x12)) + 
     &  (C12apx350pm1p*(-(kp2*1D0*(kp3 - p1p3*1D0)*
     &          ((ame**2 + p1p2)*(kp3 - kp4*1D0) + 
     &            kp2*(p1p4 - p1p3*1D0))) + 
     &       kp1**2*(-(1D0*(kp3 - p2p3*1D0)*
     &             (p2p3 - p2p4*1D0)) + 
     &          kp2*(ampi**2 - p3p4*1D0)) + 
     &       kp1*(kp2**2*(p3p4 - ampi**2*1D0) + 
     &          (ame**2 + p1p2)*(kp3 - kp4*1D0)*
     &           (kp3 - p2p3*1D0) + 
     &          kp2*(p1p3*p2p4 - p1p4*p2p3*1D0 + 
     &             kp3*(p2p3 - p1p3*1D0) + 
     &             kp4*(p1p3 - p2p3*1D0))))*2D0)/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C12a0x34x12ppp*2D0*
     &     (-(kp1**2*1D0*(p2p3 - p2p4*1D0)*
     &          (kp3 + kp4 - p2p3*1D0 - p2p4*1D0)) - 
     &       kp2*1D0*(kp3 + kp4 - p1p3*1D0 - p1p4*1D0)*
     &        ((ame**2 + p1p2)*(kp3 - kp4*1D0) + 
     &          kp2*(p1p4 - p1p3*1D0)) + 
     &       kp1*(kp3**2*(ame**2 + p1p2) - 
     &          kp4**2*(ame**2 + p1p2)*1D0 + 
     &          kp3*(-((ame**2 + p1p2)*(p2p3 + p2p4)*1D0) + 
     &             kp2*(p2p3 + p2p4 - p1p3*1D0 - p1p4*1D0)) + 
     &          kp4*((ame**2 + p1p2)*(p2p3 + p2p4) + 
     &             kp2*(p1p3 + p1p4 - p2p3*1D0 - p2p4*1D0)) + 
     &          kp2*(p1p3*p2p4 - p1p4*p2p3*1D0)*2D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C12a0x12x34ppp*2D0*
     &     (kp1**2*(p2p3 - p2p4*1D0)*
     &        (kp3 + kp4 - p2p3*1D0 - p2p4*1D0) + 
     &       kp2*(kp3 + kp4 - p1p3*1D0 - p1p4*1D0)*
     &        ((ame**2 + p1p2)*(kp3 - kp4*1D0) + 
     &          kp2*(p1p4 - p1p3*1D0)) + 
     &       kp1*(kp4**2*(ame**2 + p1p2) - 
     &          kp3**2*(ame**2 + p1p2)*1D0 + 
     &          kp4*(-((ame**2 + p1p2)*(p2p3 + p2p4)*1D0) + 
     &             kp2*(p2p3 + p2p4 - p1p3*1D0 - p1p4*1D0)) + 
     &          kp3*((ame**2 + p1p2)*(p2p3 + p2p4) + 
     &             kp2*(p1p3 + p1p4 - p2p3*1D0 - p2p4*1D0)) + 
     &          kp2*(p1p4*p2p3 - p1p3*p2p4*1D0)*2D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C00apx350pm1p*2D0*
     &     (kp1**2*(p2p3 - p2p4*1D0) + 
     &       kp1*(-((ame**2 + kp2)*1D0*
     &             (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)) + 
     &          kp4*(p1p2 + ame**2*2D0) - 
     &          kp3*1D0*(p1p2 + ame**2*2D0)) + 
     &       kp2*(ame**2*p1p4 + kp2*p1p4 + ame**2*p2p3 - 
     &          ame**2*p1p3*1D0 - kp2*p1p3*1D0 - 
     &          ame**2*p2p4*1D0 + kp3*(p1p2 + ame**2*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*2D0))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C00apx35x12pm1p*2D0*
     &     (kp1**2*(p2p3 - p2p4*1D0) + 
     &       kp1*(-((ame**2 + kp2)*1D0*
     &             (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)) + 
     &          kp4*(p1p2 + ame**2*2D0) - 
     &          kp3*1D0*(p1p2 + ame**2*2D0)) + 
     &       kp2*(ame**2*p1p4 + kp2*p1p4 + ame**2*p2p3 - 
     &          ame**2*p1p3*1D0 - kp2*p1p3*1D0 - 
     &          ame**2*p2p4*1D0 + kp3*(p1p2 + ame**2*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*2D0))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C00apx450pm1p*2D0*
     &     (kp1**2*(p2p3 - p2p4*1D0) + 
     &       kp1*(-((ame**2 + kp2)*1D0*
     &             (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)) + 
     &          kp4*(p1p2 + ame**2*2D0) - 
     &          kp3*1D0*(p1p2 + ame**2*2D0)) + 
     &       kp2*(ame**2*p1p4 + kp2*p1p4 + ame**2*p2p3 - 
     &          ame**2*p1p3*1D0 - kp2*p1p3*1D0 - 
     &          ame**2*p2p4*1D0 + kp3*(p1p2 + ame**2*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*2D0))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C00apx45x12pm1p*2D0*
     &     (kp1**2*(p2p3 - p2p4*1D0) + 
     &       kp1*(-((ame**2 + kp2)*1D0*
     &             (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)) + 
     &          kp4*(p1p2 + ame**2*2D0) - 
     &          kp3*1D0*(p1p2 + ame**2*2D0)) + 
     &       kp2*(ame**2*p1p4 + kp2*p1p4 + ame**2*p2p3 - 
     &          ame**2*p1p3*1D0 - kp2*p1p3*1D0 - 
     &          ame**2*p2p4*1D0 + kp3*(p1p2 + ame**2*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*2D0))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C1ax12x340ppp*2D0*
     &     (-(kp1**2*(ampi**2 + p3p4)*1D0*
     &          (p2p3 - p2p4*1D0)) + 
     &       kp2*(ampi**2*kp2*p1p3 - ampi**2*kp2*p1p4*1D0 + 
     &          kp2*p3p4*(p1p3 - p1p4*1D0) - 
     &          p1p3**2*p2p3*2D0 - p1p3*p1p4*p2p3*2D0 + 
     &          p1p3*p1p4*p2p4*2D0 + p1p4**2*p2p4*2D0 - 
     &          kp3*1D0*(ame**2*(ampi**2 + p3p4) + 
     &             p1p2*(ampi**2 + p3p4) - 
     &             (p1p3 + p1p4)*p2p3*2D0) + 
     &          kp4*(ame**2*(ampi**2 + p3p4) + 
     &             p1p2*(ampi**2 + p3p4) - 
     &             (p1p3 + p1p4)*p2p4*2D0)) + 
     &       kp1*(ampi**2*kp2*p1p3 + ampi**2*kp2*p2p4 - 
     &          ampi**2*kp2*p1p4*1D0 - ampi**2*kp2*p2p3*1D0 + 
     &          kp2*p3p4*
     &           (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0) + 
     &          p1p3*p2p3**2*2D0 + p1p3*p2p3*p2p4*2D0 - 
     &          p1p4*p2p3*p2p4*2D0 - p1p4*p2p4**2*2D0 + 
     &          kp3*(ame**2*(ampi**2 + p3p4) + 
     &             p1p2*(ampi**2 + p3p4) - 
     &             p1p3*(p2p3 + p2p4)*2D0) - 
     &          kp4*1D0*(ame**2*(ampi**2 + p3p4) + 
     &             p1p2*(ampi**2 + p3p4) - 
     &             p1p4*(p2p3 + p2p4)*2D0))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C1apx350pm1p*(kp1**2*
     &        (p2p3*p3p4 - ampi**2*p2p4*1D0 - 
     &          1D0*(kp3 - p2p3*1D0)*(p2p3 - p2p4*1D0) + 
     &          kp2*(ampi**2 - p3p4*1D0)) - 
     &       kp1*1D0*(ame**2*kp3*p2p3 + kp3*p1p2*p2p3 + 
     &          kp3*p1p4*p2p4 + kp3*(ame**2 + p1p2)*p3p4 - 
     &          ame**2*kp3**2*1D0 - kp3**2*p1p2*1D0 - 
     &          kp3*p1p3*p2p4*1D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0) + 
     &          kp4*(-((ame**2 + p1p2)*1D0*
     &                (ampi**2 - kp3*1D0)) + 
     &             p2p3*(p1p3 + p1p4 - ame**2*1D0 - p1p2*1D0))
     &            + kp2*(p2p3*
     &              (kp4 + p1p4 - ampi**2*1D0 - kp3*1D0) + 
     &             p3p4*(p2p4 - p1p4*1D0) + 
     &             p1p3*(ampi**2 + kp3 - kp4*1D0 - p2p4*1D0))
     &           - kp3*p1p3*p2p3*2D0 + p1p3*p2p3**2*2D0 - 
     &          p1p4*p2p3*p2p4*2D0) + 
     &       kp2*(ampi**2*kp2*p1p4 + kp2*p1p3*p1p4 - 
     &          kp3**2*(ame**2 + p1p2)*1D0 - 
     &          kp2*p1p3**2*1D0 - kp2*p1p3*p3p4*1D0 - 
     &          kp4*1D0*(ame**2*(ampi**2 + p1p3) + 
     &             p1p2*(ampi**2 + p1p3) - 
     &             p1p3*(p2p3 + p2p4)*1D0) + 
     &          p1p3**2*p2p3*2D0 - p1p3*p1p4*p2p4*2D0 + 
     &          kp3*(kp4*(ame**2 + p1p2) + 
     &             (ame**2 + p1p2)*p3p4 - 
     &             p1p4*1D0*(kp2 + p2p3 - p2p4*1D0) + 
     &             p1p3*(ame**2 + kp2 + p1p2 - p2p3*2D0)))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C1apx450pm1p*(kp1**2*
     &        (kp2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(ampi**2 + p2p4 - kp4*1D0) - 
     &          p2p4*1D0*(p2p4 + p3p4 - kp4*1D0)) + 
     &       kp1*(kp4*p1p3*p2p3 + ame**2*kp4*p2p4 + 
     &          kp4*p1p2*p2p4 + kp4*(ame**2 + p1p2)*p3p4 - 
     &          ame**2*kp4**2*1D0 - kp4**2*p1p2*1D0 - 
     &          kp4*p1p4*p2p3*1D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0) + 
     &          kp3*(-((ame**2 + p1p2)*1D0*
     &                (ampi**2 - kp4*1D0)) + 
     &             p2p4*(p1p3 + p1p4 - ame**2*1D0 - p1p2*1D0))
     &            + kp2*(p2p4*
     &              (kp3 + p1p3 - ampi**2*1D0 - kp4*1D0) + 
     &             p3p4*(p2p3 - p1p3*1D0) + 
     &             p1p4*(ampi**2 + kp4 - kp3*1D0 - p2p3*1D0))
     &           - kp4*p1p4*p2p4*2D0 - p1p3*p2p3*p2p4*2D0 + 
     &          p1p4*p2p4**2*2D0) - 
     &       kp2*1D0*(ampi**2*kp2*p1p3 - 
     &          kp4**2*(ame**2 + p1p2)*1D0 + 
     &          kp3*(kp4*(ame**2 + p1p2) - 
     &             ampi**2*(ame**2 + p1p2)*1D0 + 
     &             p1p4*(p2p3 + p2p4 - ame**2*1D0 - p1p2*1D0))
     &            - p1p4*1D0*
     &           (kp2*(p1p4 + p3p4 - p1p3*1D0) + 
     &             p1p3*p2p3*2D0 - p1p4*p2p4*2D0) + 
     &          kp4*((ame**2 + p1p2)*p3p4 + 
     &             kp2*(p1p4 - p1p3*1D0) + 
     &             p1p3*(p2p3 - p2p4*1D0) + 
     &             p1p4*(ame**2 + p1p2 - p2p4*2D0)))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C0apx45x12pm1p*(kp1**2*
     &        (kp4*p2p3 + kp2*(ampi**2 - p3p4*1D0) + 
     &          p2p4*(ampi**2 + p2p4 - kp4*1D0 - p2p3*1D0 - 
     &             p3p4*1D0)) + 
     &       kp1*(kp4**2*(ame**2 + p1p2) + 
     &          kp2**2*(p3p4 - ampi**2*1D0) + 
     &          kp2*((kp3 - kp4*1D0)*(p1p4 - p2p4*1D0) + 
     &             p2p3*(ampi**2 + p1p4 - p3p4*1D0) - 
     &             p1p3*1D0*(ampi**2 + p2p4 - p3p4*1D0)) + 
     &          p2p4*(kp3*(ame**2 + p1p2 + p1p4 - p1p3*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             (ame**2 + p1p2)*p3p4*2D0 + 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0) + 
     &          kp4*(-(kp3*(ame**2 + p1p2)*1D0) + 
     &             ame**2*(ampi**2 - p2p4*1D0 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p2p4*1D0 - p3p4*1D0) - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*2D0)))
     &        + kp2*(-(kp4**2*(ame**2 + p1p2)*1D0) + 
     &          kp4*(kp3*(ame**2 + p1p2) + 
     &             (ame**2 + p1p2)*p3p4 - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 - 
     &             p1p3*1D0*(kp2 + p2p4 - p2p3*1D0) + 
     &             p1p4*(ame**2 + kp2 + p1p2 - p2p3*2D0 + 
     &                p2p4*2D0)) - 
     &          p1p4*1D0*
     &           (kp2*(ampi**2 + p1p4 - p1p3*1D0 - 
     &                p3p4*1D0) - 
     &             1D0*(p2p3 - p2p4*1D0)*
     &              (kp3 - p1p3*2D0 + p1p4*2D0) + 
     &             ame**2*(kp3 - ampi**2*2D0 + p3p4*2D0) + 
     &             p1p2*(kp3 - ampi**2*2D0 + p3p4*2D0)))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C0apx35x12pm1p*(kp1**2*
     &        (-(p2p3**2*1D0) - kp3*p2p4*1D0 + 
     &          kp2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp3 + p2p4 + p3p4 - ampi**2*1D0)) + 
     &       kp1*(-(kp3**2*(ame**2 + p1p2)*1D0) + 
     &          kp2**2*(ampi**2 - p3p4*1D0) + 
     &          kp2*(p1p4*p2p3 - p1p3*p2p4*1D0 + 
     &             kp4*(p2p3 - p1p3*1D0) + 
     &             p3p4*(p2p4 - p1p4*1D0) + 
     &             kp3*(p1p3 - p2p3*1D0) + 
     &             ampi**2*(p1p4 - p2p4*1D0)) - 
     &          p2p3*1D0*
     &           (kp4*(ame**2 + p1p2 + p1p3 - p1p4*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             (ame**2 + p1p2)*p3p4*2D0 + 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0) + 
     &          kp3*(kp4*(ame**2 + p1p2) + 
     &             ame**2*(p2p3 + p3p4 - ampi**2*1D0) + 
     &             p1p2*(p2p3 + p3p4 - ampi**2*1D0) + 
     &             (p1p3 - p1p4*1D0)*(-(p2p4*1D0) + p2p3*2D0))
     &          ) + kp2*(kp3**2*(ame**2 + p1p2) + 
     &          kp3*(ampi**2*(ame**2 + p1p2) - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             (ame**2 + p1p2)*p3p4*1D0 + 
     &             p1p4*(kp2 + p2p3 - p2p4*1D0) - 
     &             p1p3*1D0*
     &              (ame**2 + kp2 + p1p2 + p2p3*2D0 - 
     &                p2p4*2D0)) + 
     &          p1p3*(kp2*
     &              (ampi**2 + p1p3 - p1p4*1D0 - p3p4*1D0) + 
     &             (p2p3 - p2p4*1D0)*
     &              (kp4 + p1p3*2D0 - p1p4*2D0) + 
     &             ame**2*(kp4 - ampi**2*2D0 + p3p4*2D0) + 
     &             p1p2*(kp4 - ampi**2*2D0 + p3p4*2D0)))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C12ax45x12pm1pp*2D0*
     &     (kp1**2*(-(1D0*(ampi**2 + kp3 - p2p3*1D0)*
     &             (p2p3 - p2p4*1D0)) + 
     &          kp2*(ampi**2 - p3p4*1D0)) + 
     &       kp1*(p1p4*p2p3**2 + 
     &          kp2**2*(p3p4 - ampi**2*1D0) + 
     &          kp3*(ame**2 + p1p2)*(kp3 - kp4*1D0) + 
     &          (ame**2 + p1p2)*p3p4*(kp3 + p2p4 - kp4*1D0) + 
     &          kp2*(p1p3*p2p4 - p1p4*p2p3*1D0 + 
     &             kp3*(p2p3 - p1p3*1D0) + 
     &             kp4*(p1p3 - p2p3*1D0) + 
     &             p3p4*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0))
     &           - p2p4*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*(p2p4 - kp4*1D0) + 
     &             kp3*(-(p1p4*1D0) + p1p3*2D0)) + 
     &          p2p3*(ampi**2*(ame**2 + p1p2) - 
     &             kp3*(ame**2 + p1p2 + p1p4)*1D0 - 
     &             (ame**2 + p1p2)*p3p4*1D0 + 
     &             kp4*(ame**2 + p1p2 - p1p3*1D0 + 
     &                p1p4*2D0) + p2p4*(p1p3 - p1p4*1D0)*3D0))
     &         + kp2*(-(kp3**2*(ame**2 + p1p2)*1D0) - 
     &          p1p3**2*(kp2 + p2p4)*1D0 + 
     &          (ame**2 + p1p2)*p3p4*(kp4 - p1p4*1D0) + 
     &          p1p4*(ampi**2*(ame**2 + p1p2 - kp2*1D0) + 
     &             p2p3*(p1p4 - kp4*1D0)) + 
     &          kp3*(kp4*(ame**2 + p1p2) + 
     &             p1p3*(ame**2 + kp2 + p1p2 + p2p4) - 
     &             (ame**2 + p1p2)*p3p4*1D0 - 
     &             p1p4*1D0*(kp2 + p2p4 - p2p3*2D0)) + 
     &          p1p3*(kp2*(ampi**2 + p1p4) + 
     &             (ame**2 + p1p2)*p3p4 - 
     &             ampi**2*(ame**2 + p1p2)*1D0 - 
     &             kp4*1D0*
     &              (ame**2 + p1p2 - p2p3*1D0 + p2p4*2D0) + 
     &             p1p4*(p2p4 - p2p3*1D0)*3D0))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) - 
     &  (C12ax35x12pm1pp*2D0*
     &     (kp1**2*((ampi**2 + kp4 - p2p4*1D0)*
     &           (p2p3 - p2p4*1D0) + kp2*(ampi**2 - p3p4*1D0))
     &         + kp2*(ame**2*ampi**2*p1p3 + 
     &          ampi**2*p1p2*p1p3 + ampi**2*kp2*p1p4 + 
     &          kp2*p1p3*p1p4 + p1p3**2*p2p4 - 
     &          kp4**2*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*kp2*p1p3*1D0 - 
     &          ame**2*ampi**2*p1p4*1D0 - 
     &          ampi**2*p1p2*p1p4*1D0 - kp2*p1p4**2*1D0 - 
     &          p1p4**2*p2p3*1D0 - 
     &          (ame**2 + p1p2)*p3p4*1D0*(p1p3 - p1p4*1D0) + 
     &          kp3*(kp4*(ame**2 + p1p2) + 
     &             (ame**2 + p1p2)*p3p4 - p1p3*p2p4*1D0 - 
     &             p1p4*1D0*
     &              (ame**2 + p1p2 - p2p4*1D0 + p2p3*2D0)) + 
     &          kp4*(p1p4*(ame**2 + p1p2 + p2p3) - 
     &             (ame**2 + p1p2)*p3p4*1D0 + 
     &             kp2*(p1p4 - p1p3*1D0) - 
     &             p1p3*1D0*(p2p3 - p2p4*2D0)) + 
     &          p1p3*p1p4*p2p3*3D0 - p1p3*p1p4*p2p4*3D0) + 
     &       kp1*(-(p1p4*p2p3**2*1D0) + 
     &          kp2**2*(p3p4 - ampi**2*1D0) - 
     &          kp4*(ame**2 + p1p2)*1D0*(kp3 - kp4*1D0) - 
     &          (ame**2 + p1p2)*p3p4*1D0*
     &           (kp3 + p2p4 - kp4*1D0) + 
     &          kp2*(p1p4*p2p3 - p1p3*(p2p4 + p3p4)*1D0 + 
     &             kp4*(p2p4 - p1p4*1D0) + 
     &             kp3*(p1p4 - p2p4*1D0) + 
     &             p3p4*(p1p4 + p2p3 - p2p4*1D0)) + 
     &          p2p4*(ampi**2*(ame**2 + p1p2) + p1p3*p2p4 - 
     &             kp4*(ame**2 + p1p2 + p1p3)*1D0 + 
     &             kp3*(ame**2 + p1p2 - p1p4*1D0 + p1p3*2D0))
     &           + p2p3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             kp4*(p1p3 - p1p4*2D0) - p1p3*p2p4*3D0 + 
     &             p1p4*(kp3 + p2p4*3D0)))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C00a0x12x34ppp*(kp1**2*(p2p4 - p2p3*1D0) + 
     &       kp2*(ame**2*p1p3 + kp2*p1p3 + ame**2*p2p4 - 
     &          ame**2*p1p4*1D0 - kp2*p1p4*1D0 - 
     &          ame**2*p2p3*1D0 + kp4*(p1p2 + ame**2*2D0) - 
     &          kp3*1D0*(p1p2 + ame**2*2D0)) + 
     &       kp1*((ame**2 + kp2)*
     &           (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0) + 
     &          kp3*(p1p2 + ame**2*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*2D0)))*4D0)/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D00ap0x12px35x34m1ppp*2D0*
     &     (kp1**2*(p2p4 - p2p3*1D0) + 
     &       kp2*(ame**2*p1p3 + kp2*p1p3 + ame**2*p2p4 - 
     &          ame**2*p1p4*1D0 - kp2*p1p4*1D0 - 
     &          ame**2*p2p3*1D0 + kp4*(p1p2 + ame**2*2D0) - 
     &          kp3*1D0*(p1p2 + ame**2*2D0)) + 
     &       kp1*((ame**2 + kp2)*
     &           (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0) + 
     &          kp3*(p1p2 + ame**2*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*2D0)))*
     &     (m1**2 + p3p4*4D0))/(kp1*kp2*(ampi**2 + p3p4)*x12)
     &   + (D00ap0x12px45x34m1ppp*2D0*
     &     (kp1**2*(p2p4 - p2p3*1D0) + 
     &       kp2*(ame**2*p1p3 + kp2*p1p3 + ame**2*p2p4 - 
     &          ame**2*p1p4*1D0 - kp2*p1p4*1D0 - 
     &          ame**2*p2p3*1D0 + kp4*(p1p2 + ame**2*2D0) - 
     &          kp3*1D0*(p1p2 + ame**2*2D0)) + 
     &       kp1*((ame**2 + kp2)*
     &           (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0) + 
     &          kp3*(p1p2 + ame**2*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*2D0)))*
     &     (m1**2 + p3p4*4D0))/(kp1*kp2*(ampi**2 + p3p4)*x12)
     &   + (D22ap0x12px45x34m1ppp*
     &     (-(kp1**2*1D0*
     &          (p2p3*x45 + kp2*(ampi**2 - p3p4*1D0)*2D0 + 
     &            p2p4*(ampi**2 - kp4*2D0 - p2p3*2D0 + 
     &               p2p4*2D0 - p3p4*2D0))) + 
     &       kp1*(kp3*(ame**2 + p1p2)*x45 - 
     &          kp4**2*(ame**2 + p1p2)*2D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0)*2D0 + 
     &          kp4*(-(ampi**2*(ame**2 + p1p2)*1D0) + 
     &             p2p4*2D0*
     &              (ame**2 + p1p2 - p1p3*1D0 + p1p4*2D0)) - 
     &          p2p4*2D0*
     &           (kp3*(ame**2 + p1p2 + p1p4) + 
     &             (ame**2 + p1p2)*p3p4 - 
     &             ampi**2*(ame**2 + p1p2)*1D0 - 
     &             p1p3*p2p4*1D0 - p1p4*1D0*(p2p3 - p2p4*2D0))
     &            + kp2*(x35*(p2p4 - p1p4*1D0) + 
     &             kp4*(p1p4 - p2p4*1D0)*2D0 - 
     &             p2p3*1D0*(ampi**2 + p1p4*2D0) + 
     &             p1p3*(ampi**2 + p2p4*2D0))) + 
     &       kp2*(kp4**2*x12 - kp3*(ame**2 + p1p2)*x45*1D0 + 
     &          p1p4*2D0*
     &           (kp3*(ame**2 + p1p2 + p2p4) + 
     &             (ame**2 + p1p2)*p3p4 - 
     &             ampi**2*(ame**2 + p1p2)*1D0 - 
     &             p1p4*p2p3*1D0 - p2p4*1D0*(p1p3 - p1p4*2D0))
     &            + kp4*(ame**2*ampi**2 + 
     &             kp2*(p1p3 - p1p4*1D0)*2D0 + 
     &             p1p2*(ampi**2 - p1p4*2D0) - 
     &             p1p4*2D0*(ame**2 - p2p3*1D0 + p2p4*2D0)) + 
     &          kp2*(p1p3*(ampi**2 - p1p4*2D0) + 
     &             p1p4*(ampi**2 + p1p4*2D0 - p3p4*2D0))))*
     &     (m1**2 + p3p4*4D0))/(kp1*kp2*(ampi**2 + p3p4)*x12)
     &   + (D0ap0x12px45x34m1ppp*
     &     (kp1**2*(-(kp4*p2p3*1D0) + 
     &          kp2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp4 + p2p3 + p3p4 - ampi**2*1D0 - 
     &             p2p4*1D0)) + 
     &       kp1*(-(kp4**2*(ame**2 + p1p2)*1D0) + 
     &          kp2**2*(ampi**2 - p3p4*1D0) + 
     &          kp2*(-(1D0*(kp3 - kp4*1D0)*
     &                (p1p4 - p2p4*1D0)) - 
     &             p2p3*1D0*(ampi**2 + p1p4 - p3p4*1D0) + 
     &             p1p3*(ampi**2 + p2p4 - p3p4*1D0)) + 
     &          p2p4*(ampi**2*x12 - 
     &             kp3*1D0*
     &              (ame**2 + p1p2 + p1p4 - p1p3*1D0) - 
     &             (ame**2 + p1p2)*p3p4*2D0 - 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0) + 
     &          kp4*(kp3*(ame**2 + p1p2) + 
     &             ame**2*(p2p4 + p3p4 - ampi**2*1D0) + 
     &             p1p2*(p2p4 + p3p4 - ampi**2*1D0) + 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*2D0))) + 
     &       kp2*(kp4**2*(ame**2 + p1p2) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp3*(ame**2 + p1p2)*1D0 - 
     &             (ame**2 + p1p2)*p3p4*1D0 + 
     &             p1p3*(kp2 + p2p4 - p2p3*1D0) - 
     &             p1p4*1D0*
     &              (ame**2 + kp2 + p1p2 - p2p3*2D0 + 
     &                p2p4*2D0)) + 
     &          p1p4*(kp2*
     &              (ampi**2 + p1p4 - p1p3*1D0 - p3p4*1D0) - 
     &             1D0*(p2p3 - p2p4*1D0)*
     &              (kp3 - p1p3*2D0 + p1p4*2D0) + 
     &             ame**2*(kp3 - ampi**2*2D0 + p3p4*2D0) + 
     &             p1p2*(kp3 - ampi**2*2D0 + p3p4*2D0))))*
     &     (m1**2 + p3p4*4D0))/(kp1*kp2*(ampi**2 + p3p4)*x12)
     &   + (D0ap0x12px35x34m1ppp*
     &     (kp1**2*(p2p3**2 + kp3*p2p4 + 
     &          kp2*(ampi**2 - p3p4*1D0) + 
     &          p2p3*(ampi**2 - kp3*1D0 - p2p4*1D0 - p3p4*1D0)
     &          ) + kp1*(kp3**2*(ame**2 + p1p2) + 
     &          kp2**2*(p3p4 - ampi**2*1D0) + 
     &          kp2*(p1p3*p2p4 - p1p4*p2p3*1D0 + 
     &             kp3*(p2p3 - p1p3*1D0) + 
     &             ampi**2*(p2p4 - p1p4*1D0) + 
     &             kp4*(p1p3 - p2p3*1D0) + 
     &             p3p4*(p1p4 - p2p4*1D0)) + 
     &          p2p3*(p3p4*x12 + 
     &             kp4*(ame**2 + p1p2 + p1p3 - p1p4*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0) + 
     &          kp3*(ampi**2*(ame**2 + p1p2) - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             (ame**2 + p1p2)*p3p4*1D0 + 
     &             p2p4*(p1p3 - p1p4*1D0) - 
     &             p2p3*1D0*
     &              (ame**2 + p1p2 + p1p3*2D0 - p1p4*2D0))) + 
     &       kp2*(-(kp3**2*(ame**2 + p1p2)*1D0) + 
     &          kp3*(kp4*(ame**2 + p1p2) + 
     &             (ame**2 + p1p2)*p3p4 - 
     &             ampi**2*(ame**2 + p1p2)*1D0 - 
     &             p1p4*1D0*(kp2 + p2p3 - p2p4*1D0) + 
     &             p1p3*(ame**2 + kp2 + p1p2 + p2p3*2D0 - 
     &                p2p4*2D0)) - 
     &          p1p3*1D0*
     &           (kp2*(ampi**2 + p1p3 - p1p4*1D0 - 
     &                p3p4*1D0) + 
     &             (p2p3 - p2p4*1D0)*
     &              (kp4 + p1p3*2D0 - p1p4*2D0) + 
     &             ame**2*(kp4 - ampi**2*2D0 + p3p4*2D0) + 
     &             p1p2*(kp4 - ampi**2*2D0 + p3p4*2D0))))*
     &     (m1**2 + p3p4*4D0))/(kp1*kp2*(ampi**2 + p3p4)*x12)
     &   + (D13ap0x12px35x34m1ppp*2D0*
     &     (ampi**2*kp2**2*(p1p4 - p1p3*1D0) + 
     &       kp1*(-(p1p4*p2p3**2*1D0) - 
     &          (ame**2 + p1p2)*p3p4*1D0*
     &           (kp3 + p2p4 - kp4*1D0) + 
     &          p2p4*(ame**2*ampi**2 + 
     &             ampi**2*(p1p2 - kp1*1D0) + 
     &             p1p3*(p2p4 - kp4*1D0) + 
     &             kp3*(-(p1p4*1D0) + p1p3*2D0)) + 
     &          p2p3*(ampi**2*kp1 + kp3*p1p4 + 
     &             ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             kp4*(p1p3 - p1p4*2D0) + 
     &             p2p4*(p1p4 - p1p3*1D0)*3D0)) + 
     &       kp2*(p1p3**2*p2p4 - p1p4**2*p2p3*1D0 + 
     &          p3p4*((ame**2 + p1p2)*(kp3 - kp4*1D0) + 
     &             kp1*(p2p3 - p2p4*1D0)) + 
     &          p1p4*(kp4*p2p3 + (ame**2 + kp1)*p3p4 - 
     &             ame**2*ampi**2*1D0 + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             kp3*(p2p4 - p2p3*2D0)) + 
     &          p1p3*(ame**2*ampi**2 - kp3*p2p4*1D0 - 
     &             (ame**2 + kp1)*p3p4*1D0 + 
     &             p1p2*(ampi**2 - p3p4*1D0) - 
     &             kp4*1D0*(p2p3 - p2p4*2D0) + 
     &             p1p4*(p2p3 - p2p4*1D0)*3D0)))*
     &     (m1**2 + p3p4*4D0))/(kp1*kp2*(ampi**2 + p3p4)*x12)
     &   + (D13ap0x12px45x34m1ppp*2D0*
     &     (ampi**2*kp2**2*(p1p4 - p1p3*1D0) + 
     &       kp1*(-(p1p4*p2p3**2*1D0) - 
     &          (ame**2 + p1p2)*p3p4*1D0*
     &           (kp3 + p2p4 - kp4*1D0) + 
     &          p2p4*(ame**2*ampi**2 + 
     &             ampi**2*(p1p2 - kp1*1D0) + 
     &             p1p3*(p2p4 - kp4*1D0) + 
     &             kp3*(-(p1p4*1D0) + p1p3*2D0)) + 
     &          p2p3*(ampi**2*kp1 + kp3*p1p4 + 
     &             ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             kp4*(p1p3 - p1p4*2D0) + 
     &             p2p4*(p1p4 - p1p3*1D0)*3D0)) + 
     &       kp2*(p1p3**2*p2p4 - p1p4**2*p2p3*1D0 + 
     &          p3p4*((ame**2 + p1p2)*(kp3 - kp4*1D0) + 
     &             kp1*(p2p3 - p2p4*1D0)) + 
     &          p1p4*(kp4*p2p3 + (ame**2 + kp1)*p3p4 - 
     &             ame**2*ampi**2*1D0 + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             kp3*(p2p4 - p2p3*2D0)) + 
     &          p1p3*(ame**2*ampi**2 - kp3*p2p4*1D0 - 
     &             (ame**2 + kp1)*p3p4*1D0 + 
     &             p1p2*(ampi**2 - p3p4*1D0) - 
     &             kp4*1D0*(p2p3 - p2p4*2D0) + 
     &             p1p4*(p2p3 - p2p4*1D0)*3D0)))*
     &     (m1**2 + p3p4*4D0))/(kp1*kp2*(ampi**2 + p3p4)*x12)
     &   + (D23ap0x12px45x34m1ppp*2D0*
     &     (kp1**2*(kp2*(p3p4 - ampi**2*1D0) + 
     &          (ampi**2 + kp3 - p2p3*1D0)*(p2p3 - p2p4*1D0))
     &        + kp1*(-(p1p4*p2p3**2*1D0) - 
     &          kp3*(ame**2 + p1p2)*1D0*(kp3 - kp4*1D0) - 
     &          (ame**2 + p1p2)*p3p4*1D0*
     &           (kp3 + p2p4 - kp4*1D0) + 
     &          kp2**2*(ampi**2 - p3p4*1D0) + 
     &          kp2*(p1p4*p2p3 - p1p3*(p2p4 + p3p4)*1D0 + 
     &             kp4*(p2p3 - p1p3*1D0) + 
     &             kp3*(p1p3 - p2p3*1D0) + 
     &             p3p4*(p1p4 + p2p3 - p2p4*1D0)) + 
     &          p2p4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*(p2p4 - kp4*1D0) + 
     &             kp3*(-(p1p4*1D0) + p1p3*2D0)) + 
     &          p2p3*(kp3*(ame**2 + p1p2 + p1p4) + 
     &             (ame**2 + p1p2)*p3p4 - 
     &             ampi**2*(ame**2 + p1p2)*1D0 - 
     &             kp4*1D0*
     &              (ame**2 + p1p2 - p1p3*1D0 + p1p4*2D0) + 
     &             p2p4*(p1p4 - p1p3*1D0)*3D0)) + 
     &       kp2*(kp3**2*(ame**2 + p1p2) + 
     &          p1p3**2*(kp2 + p2p4) + 
     &          kp3*(ame**2 + p1p2)*p3p4 + 
     &          (ame**2 + p1p2)*p3p4*(p1p4 - kp4*1D0) + 
     &          p1p4*(-(ampi**2*1D0*
     &                (ame**2 + p1p2 - kp2*1D0)) + 
     &             p2p3*(kp4 - p1p4*1D0)) - 
     &          kp3*1D0*(kp4*(ame**2 + p1p2) + 
     &             p1p3*(ame**2 + kp2 + p1p2 + p2p4) - 
     &             p1p4*1D0*(kp2 + p2p4 - p2p3*2D0)) + 
     &          p1p3*(ampi**2*(ame**2 + p1p2) - 
     &             kp2*(ampi**2 + p1p4)*1D0 - 
     &             (ame**2 + p1p2)*p3p4*1D0 + 
     &             kp4*(ame**2 + p1p2 - p2p3*1D0 + 
     &                p2p4*2D0) + p1p4*(p2p3 - p2p4*1D0)*3D0))
     &       )*(m1**2 + p3p4*4D0))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D23ap0x12px35x34m1ppp*2D0*
     &     (kp1**2*((ampi**2 + kp4 - p2p4*1D0)*
     &           (p2p3 - p2p4*1D0) + kp2*(ampi**2 - p3p4*1D0))
     &         + kp2*(ame**2*ampi**2*p1p3 + 
     &          ampi**2*p1p2*p1p3 + ampi**2*kp2*p1p4 + 
     &          kp2*p1p3*p1p4 + p1p3**2*p2p4 - 
     &          kp4**2*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*kp2*p1p3*1D0 - 
     &          ame**2*ampi**2*p1p4*1D0 - 
     &          ampi**2*p1p2*p1p4*1D0 - kp2*p1p4**2*1D0 - 
     &          p1p4**2*p2p3*1D0 - 
     &          (ame**2 + p1p2)*p3p4*1D0*(p1p3 - p1p4*1D0) + 
     &          kp3*(kp4*(ame**2 + p1p2) + 
     &             (ame**2 + p1p2)*p3p4 - p1p3*p2p4*1D0 - 
     &             p1p4*1D0*
     &              (ame**2 + p1p2 - p2p4*1D0 + p2p3*2D0)) + 
     &          kp4*(p1p4*(ame**2 + p1p2 + p2p3) - 
     &             (ame**2 + p1p2)*p3p4*1D0 + 
     &             kp2*(p1p4 - p1p3*1D0) - 
     &             p1p3*1D0*(p2p3 - p2p4*2D0)) + 
     &          p1p3*p1p4*p2p3*3D0 - p1p3*p1p4*p2p4*3D0) + 
     &       kp1*(-(p1p4*p2p3**2*1D0) + 
     &          kp2**2*(p3p4 - ampi**2*1D0) - 
     &          kp4*(ame**2 + p1p2)*1D0*(kp3 - kp4*1D0) - 
     &          (ame**2 + p1p2)*p3p4*1D0*
     &           (kp3 + p2p4 - kp4*1D0) + 
     &          kp2*(p1p4*p2p3 - p1p3*(p2p4 + p3p4)*1D0 + 
     &             kp4*(p2p4 - p1p4*1D0) + 
     &             kp3*(p1p4 - p2p4*1D0) + 
     &             p3p4*(p1p4 + p2p3 - p2p4*1D0)) + 
     &          p2p4*(ampi**2*(ame**2 + p1p2) + p1p3*p2p4 - 
     &             kp4*(ame**2 + p1p2 + p1p3)*1D0 + 
     &             kp3*(ame**2 + p1p2 - p1p4*1D0 + p1p3*2D0))
     &           + p2p3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             kp4*(p1p3 - p1p4*2D0) - p1p3*p2p4*3D0 + 
     &             p1p4*(kp3 + p2p4*3D0))))*(m1**2 + p3p4*4D0)
     &     )/(kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C11apx350pm1p*(-(kp1**2*1D0*
     &          (ampi**2*p2p4 + p2p3*(ampi**2 - p3p4*2D0))) + 
     &       kp2*(ampi**2*kp2*p1p3 + ampi**2*kp2*p1p4 - 
     &          ame**2*ampi**2*p1p3*2D0 - 
     &          ampi**2*p1p2*p1p3*2D0 - p1p3*p1p4*p2p3*2D0 - 
     &          p1p3**2*p2p4*2D0 + 
     &          p1p3*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 - 
     &          kp4*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*2D0) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*2D0*(p2p4 - p2p3*2D0)) + 
     &          p1p3**2*p2p3*4D0) + 
     &       kp1*(ampi**2*kp2*p1p4 + ampi**2*kp2*p2p3 - 
     &          ampi**2*kp2*p1p3*1D0 - ampi**2*kp2*p2p4*1D0 + 
     &          ame**2*ampi**2*p2p3*2D0 + 
     &          ampi**2*p1p2*p2p3*2D0 + p1p4*p2p3**2*2D0 + 
     &          p1p3*p2p3*p2p4*2D0 - 
     &          (ame**2 + p1p2)*p2p3*p3p4*2D0 + 
     &          kp4*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*2D0) + 
     &          kp3*(-(ampi**2*(ame**2 + p1p2)*1D0) + 
     &             p2p3*2D0*(-(p1p4*1D0) + p1p3*2D0)) - 
     &          p1p3*p2p3**2*4D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C22ax45px12pm1p*
     &     (-(kp1**2*1D0*
     &          (ampi**2*p2p4 + p2p3*(ampi**2 - p3p4*2D0))) + 
     &       kp2*(ampi**2*kp2*p1p3 + ampi**2*kp2*p1p4 - 
     &          ame**2*ampi**2*p1p3*2D0 - 
     &          ampi**2*p1p2*p1p3*2D0 - p1p3*p1p4*p2p3*2D0 - 
     &          p1p3**2*p2p4*2D0 + 
     &          p1p3*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 - 
     &          kp4*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*2D0) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*2D0*(p2p4 - p2p3*2D0)) + 
     &          p1p3**2*p2p3*4D0) + 
     &       kp1*(ampi**2*kp2*p1p4 + ampi**2*kp2*p2p3 - 
     &          ampi**2*kp2*p1p3*1D0 - ampi**2*kp2*p2p4*1D0 + 
     &          ame**2*ampi**2*p2p3*2D0 + 
     &          ampi**2*p1p2*p2p3*2D0 + p1p4*p2p3**2*2D0 + 
     &          p1p3*p2p3*p2p4*2D0 - 
     &          (ame**2 + p1p2)*p2p3*p3p4*2D0 + 
     &          kp4*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*2D0) + 
     &          kp3*(-(ampi**2*(ame**2 + p1p2)*1D0) + 
     &             p2p3*2D0*(-(p1p4*1D0) + p1p3*2D0)) - 
     &          p1p3*p2p3**2*4D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C22apx35x12pm1p*
     &     (-(kp1**2*1D0*
     &          (p2p3**2*2D0 + kp2*(ampi**2 - p3p4*1D0)*2D0 + 
     &            p2p4*(ampi**2 + kp3*2D0) + 
     &            p2p3*(ampi**2 - kp3*2D0 - p2p4*2D0 - 
     &               p3p4*2D0))) + 
     &       kp2*(ampi**2*kp2*p1p3 + ampi**2*kp2*p1p4 + 
     &          kp3**2*(ame**2 + p1p2)*2D0 - 
     &          ame**2*ampi**2*p1p3*2D0 - 
     &          ampi**2*p1p2*p1p3*2D0 + kp2*p1p3**2*2D0 - 
     &          kp2*p1p3*p1p4*2D0 - p1p3*p1p4*p2p3*2D0 - 
     &          p1p3**2*p2p4*2D0 + 
     &          p1p3*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 + 
     &          kp4*(p1p3*p2p3*2D0 - 
     &             (ame**2 + p1p2)*1D0*(ampi**2 - p1p3*2D0))
     &           + kp3*(ampi**2*(ame**2 + p1p2) - 
     &             kp4*(ame**2 + p1p2)*2D0 + kp2*p1p4*2D0 - 
     &             p1p3*2D0*
     &              (ame**2 + kp2 + p1p2 - p2p4*1D0 + 
     &                p2p3*2D0)) + p1p3**2*p2p3*4D0) + 
     &       kp1*(ame**2*ampi**2*kp4 + ampi**2*kp4*p1p2 - 
     &          kp3**2*(ame**2 + p1p2)*2D0 + 
     &          ame**2*ampi**2*p2p3*2D0 - 
     &          ame**2*kp4*p2p3*2D0 + ampi**2*p1p2*p2p3*2D0 - 
     &          kp4*p1p2*p2p3*2D0 - kp4*p1p3*p2p3*2D0 + 
     &          p1p4*p2p3**2*2D0 + p1p3*p2p3*p2p4*2D0 - 
     &          (ame**2 + p1p2)*p2p3*p3p4*2D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0)*2D0 + 
     &          kp3*(-((ame**2 + p1p2)*1D0*
     &                (ampi**2 - kp4*2D0)) + 
     &             p2p3*2D0*
     &              (ame**2 + p1p2 - p1p4*1D0 + p1p3*2D0)) + 
     &          kp2*(-(ampi**2*p2p4*1D0) + 
     &             p2p3*(ampi**2 - kp3*2D0 + kp4*2D0) + 
     &             p1p4*(ampi**2 + p2p3*2D0) - 
     &             p1p3*1D0*
     &              (ampi**2 - kp3*2D0 + kp4*2D0 + p2p4*2D0))
     &           - p1p3*p2p3**2*4D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D11ap0x12px35x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (kp1**2*(ampi**2*p2p4 + 
     &          p2p3*(ampi**2 - p3p4*2D0)) + 
     &       kp2*(-(ampi**2*kp2*p1p3*1D0) - 
     &          ampi**2*kp2*p1p4*1D0 + 
     &          ame**2*ampi**2*p1p3*2D0 + 
     &          ampi**2*p1p2*p1p3*2D0 + p1p3*p1p4*p2p3*2D0 + 
     &          p1p3**2*p2p4*2D0 - 
     &          p1p3*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 + 
     &          kp4*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*2D0) - 
     &          kp3*1D0*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*2D0*(p2p4 - p2p3*2D0)) - 
     &          p1p3**2*p2p3*4D0) + 
     &       kp1*(ampi**2*kp2*p1p3 + ampi**2*kp2*p2p4 + 
     &          p2p3*p3p4*x12 - ampi**2*kp2*p1p4*1D0 - 
     &          ampi**2*kp2*p2p3*1D0 - 
     &          ame**2*ampi**2*p2p3*2D0 - 
     &          ampi**2*p1p2*p2p3*2D0 - p1p4*p2p3**2*2D0 - 
     &          p1p3*p2p3*p2p4*2D0 - 
     &          kp4*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*2D0) + 
     &          kp3*(ampi**2*(ame**2 + p1p2) + 
     &             p2p3*2D0*(p1p4 - p1p3*2D0)) + 
     &          p1p3*p2p3**2*4D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D33ap0x12px45x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (kp1**2*(ampi**2*p2p4 + 
     &          p2p3*(ampi**2 - p3p4*2D0)) + 
     &       kp2*(-(ampi**2*kp2*p1p3*1D0) - 
     &          ampi**2*kp2*p1p4*1D0 + 
     &          ame**2*ampi**2*p1p3*2D0 + 
     &          ampi**2*p1p2*p1p3*2D0 + p1p3*p1p4*p2p3*2D0 + 
     &          p1p3**2*p2p4*2D0 - 
     &          p1p3*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 + 
     &          kp4*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*2D0) - 
     &          kp3*1D0*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*2D0*(p2p4 - p2p3*2D0)) - 
     &          p1p3**2*p2p3*4D0) + 
     &       kp1*(ampi**2*kp2*p1p3 + ampi**2*kp2*p2p4 + 
     &          p2p3*p3p4*x12 - ampi**2*kp2*p1p4*1D0 - 
     &          ampi**2*kp2*p2p3*1D0 - 
     &          ame**2*ampi**2*p2p3*2D0 - 
     &          ampi**2*p1p2*p2p3*2D0 - p1p4*p2p3**2*2D0 - 
     &          p1p3*p2p3*p2p4*2D0 - 
     &          kp4*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*2D0) + 
     &          kp3*(ampi**2*(ame**2 + p1p2) + 
     &             p2p3*2D0*(p1p4 - p1p3*2D0)) + 
     &          p1p3*p2p3**2*4D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D22ap0x12px35x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (kp1**2*(p2p4*x35 + kp2*(ampi**2 - p3p4*1D0)*2D0 + 
     &          p2p3*(ampi**2 - kp3*2D0 + p2p3*2D0 - 
     &             p2p4*2D0 - p3p4*2D0)) - 
     &       kp2*1D0*(ampi**2*kp2*p1p3 + ampi**2*kp2*p1p4 + 
     &          kp3**2*x12 - ame**2*ampi**2*p1p3*2D0 - 
     &          ampi**2*p1p2*p1p3*2D0 + kp2*p1p3**2*2D0 - 
     &          kp2*p1p3*p1p4*2D0 - p1p3*p1p4*p2p3*2D0 - 
     &          p1p3**2*p2p4*2D0 + 
     &          p1p3*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 + 
     &          kp4*(p1p3*p2p3*2D0 - 
     &             (ame**2 + p1p2)*1D0*(ampi**2 - p1p3*2D0))
     &           + kp3*(ampi**2*(ame**2 + p1p2) - 
     &             kp4*(ame**2 + p1p2)*2D0 + kp2*p1p4*2D0 - 
     &             p1p3*2D0*
     &              (ame**2 + kp2 + p1p2 - p2p4*1D0 + 
     &                p2p3*2D0)) + p1p3**2*p2p3*4D0) + 
     &       kp1*(kp3**2*x12 + p2p3*p3p4*x12 - 
     &          ame**2*ampi**2*kp4*1D0 - 
     &          ampi**2*kp4*p1p2*1D0 - 
     &          ame**2*ampi**2*p2p3*2D0 + 
     &          ame**2*kp4*p2p3*2D0 - ampi**2*p1p2*p2p3*2D0 + 
     &          kp4*p1p2*p2p3*2D0 + kp4*p1p3*p2p3*2D0 - 
     &          p1p4*p2p3**2*2D0 - p1p3*p2p3*p2p4*2D0 + 
     &          kp2**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &          kp3*((ame**2 + p1p2)*(ampi**2 - kp4*2D0) - 
     &             p2p3*2D0*
     &              (ame**2 + p1p2 - p1p4*1D0 + p1p3*2D0)) + 
     &          kp2*(ampi**2*(p2p4 - p1p4*1D0) + 
     &             kp3*(p2p3 - p1p3*1D0)*2D0 - 
     &             p2p3*1D0*(x45 + p1p4*2D0) + 
     &             p1p3*(x45 + p2p4*2D0)) + p1p3*p2p3**2*4D0))
     &     )/(kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D11ap0x12px45x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (-(kp1**2*1D0*
     &          (ampi**2*p2p3 + p2p4*(ampi**2 - p3p4*2D0))) + 
     &       kp2*(ampi**2*kp2*p1p3 + ampi**2*kp2*p1p4 - 
     &          ame**2*ampi**2*p1p4*2D0 - 
     &          ampi**2*p1p2*p1p4*2D0 - p1p4**2*p2p3*2D0 - 
     &          p1p3*p1p4*p2p4*2D0 + 
     &          p1p4*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p4*p2p4*2D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p4*2D0*(p2p3 - p2p4*2D0)) + 
     &          p1p4**2*p2p4*4D0) + 
     &       kp1*(ampi**2*kp2*p1p3 + ampi**2*kp2*p2p4 - 
     &          ampi**2*kp2*p1p4*1D0 - ampi**2*kp2*p2p3*1D0 + 
     &          ame**2*ampi**2*p2p4*2D0 + 
     &          ampi**2*p1p2*p2p4*2D0 + p1p4*p2p3*p2p4*2D0 + 
     &          p1p3*p2p4**2*2D0 - 
     &          (ame**2 + p1p2)*p2p4*p3p4*2D0 + 
     &          kp3*(ampi**2*(ame**2 + p1p2) - 
     &             p1p4*p2p4*2D0) + 
     &          kp4*(-(ampi**2*(ame**2 + p1p2)*1D0) - 
     &             p2p4*2D0*(p1p3 - p1p4*2D0)) - 
     &          p1p4*p2p4**2*4D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D33ap0x12px35x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (-(kp1**2*1D0*
     &          (ampi**2*p2p3 + p2p4*(ampi**2 - p3p4*2D0))) + 
     &       kp2*(ampi**2*kp2*p1p3 + ampi**2*kp2*p1p4 - 
     &          ame**2*ampi**2*p1p4*2D0 - 
     &          ampi**2*p1p2*p1p4*2D0 - p1p4**2*p2p3*2D0 - 
     &          p1p3*p1p4*p2p4*2D0 + 
     &          p1p4*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p4*p2p4*2D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p4*2D0*(p2p3 - p2p4*2D0)) + 
     &          p1p4**2*p2p4*4D0) + 
     &       kp1*(ampi**2*kp2*p1p3 + ampi**2*kp2*p2p4 - 
     &          ampi**2*kp2*p1p4*1D0 - ampi**2*kp2*p2p3*1D0 + 
     &          ame**2*ampi**2*p2p4*2D0 + 
     &          ampi**2*p1p2*p2p4*2D0 + p1p4*p2p3*p2p4*2D0 + 
     &          p1p3*p2p4**2*2D0 - 
     &          (ame**2 + p1p2)*p2p4*p3p4*2D0 + 
     &          kp3*(ampi**2*(ame**2 + p1p2) - 
     &             p1p4*p2p4*2D0) + 
     &          kp4*(-(ampi**2*(ame**2 + p1p2)*1D0) - 
     &             p2p4*2D0*(p1p3 - p1p4*2D0)) - 
     &          p1p4*p2p4**2*4D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C11apx450pm1p*(kp1**2*
     &        (ampi**2*p2p3 + p2p4*(ampi**2 - p3p4*2D0)) + 
     &       kp2*(-(ampi**2*kp2*p1p3*1D0) - 
     &          ampi**2*kp2*p1p4*1D0 + 
     &          ame**2*ampi**2*p1p4*2D0 + 
     &          ampi**2*p1p2*p1p4*2D0 + p1p4**2*p2p3*2D0 + 
     &          p1p3*p1p4*p2p4*2D0 - 
     &          p1p4*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 + 
     &          kp3*(ampi**2*(ame**2 + p1p2) - 
     &             p1p4*p2p4*2D0) - 
     &          kp4*1D0*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p4*2D0*(p2p3 - p2p4*2D0)) - 
     &          p1p4**2*p2p4*4D0) + 
     &       kp1*(ampi**2*kp2*p1p4 + ampi**2*kp2*p2p3 - 
     &          ampi**2*kp2*p1p3*1D0 - ampi**2*kp2*p2p4*1D0 - 
     &          ame**2*ampi**2*p2p4*2D0 - 
     &          ampi**2*p1p2*p2p4*2D0 - p1p4*p2p3*p2p4*2D0 - 
     &          p1p3*p2p4**2*2D0 + 
     &          (ame**2 + p1p2)*p2p4*p3p4*2D0 - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p4*p2p4*2D0) + 
     &          kp4*(ampi**2*(ame**2 + p1p2) + 
     &             p2p4*2D0*(p1p3 - p1p4*2D0)) + 
     &          p1p4*p2p4**2*4D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C22ax35px12pm1p*
     &     (kp1**2*(ampi**2*p2p3 + 
     &          p2p4*(ampi**2 - p3p4*2D0)) + 
     &       kp2*(-(ampi**2*kp2*p1p3*1D0) - 
     &          ampi**2*kp2*p1p4*1D0 + 
     &          ame**2*ampi**2*p1p4*2D0 + 
     &          ampi**2*p1p2*p1p4*2D0 + p1p4**2*p2p3*2D0 + 
     &          p1p3*p1p4*p2p4*2D0 - 
     &          p1p4*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 + 
     &          kp3*(ampi**2*(ame**2 + p1p2) - 
     &             p1p4*p2p4*2D0) - 
     &          kp4*1D0*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p4*2D0*(p2p3 - p2p4*2D0)) - 
     &          p1p4**2*p2p4*4D0) + 
     &       kp1*(ampi**2*kp2*p1p4 + ampi**2*kp2*p2p3 - 
     &          ampi**2*kp2*p1p3*1D0 - ampi**2*kp2*p2p4*1D0 - 
     &          ame**2*ampi**2*p2p4*2D0 - 
     &          ampi**2*p1p2*p2p4*2D0 - p1p4*p2p3*p2p4*2D0 - 
     &          p1p3*p2p4**2*2D0 + 
     &          (ame**2 + p1p2)*p2p4*p3p4*2D0 - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p4*p2p4*2D0) + 
     &          kp4*(ampi**2*(ame**2 + p1p2) + 
     &             p2p4*2D0*(p1p3 - p1p4*2D0)) + 
     &          p1p4*p2p4**2*4D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C22apx45x12pm1p*
     &     (kp1**2*(kp2*(ampi**2 - p3p4*1D0)*2D0 + 
     &          p2p3*(ampi**2 + kp4*2D0 - p2p4*2D0) + 
     &          p2p4*(ampi**2 - kp4*2D0 + p2p4*2D0 - p3p4*2D0)
     &          ) - kp2*1D0*
     &        (ampi**2*kp2*p1p3 + ampi**2*kp2*p1p4 + 
     &          kp4**2*(ame**2 + p1p2)*2D0 - 
     &          ame**2*ampi**2*p1p4*2D0 - 
     &          ampi**2*p1p2*p1p4*2D0 - kp2*p1p3*p1p4*2D0 + 
     &          kp2*p1p4**2*2D0 - p1p4**2*p2p3*2D0 - 
     &          p1p3*p1p4*p2p4*2D0 + 
     &          p1p4*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 + 
     &          kp3*(p1p4*p2p4*2D0 - 
     &             (ame**2 + p1p2)*1D0*
     &              (ampi**2 + kp4*2D0 - p1p4*2D0)) + 
     &          kp4*(ame**2*ampi**2 + 
     &             kp2*(p1p3 - p1p4*1D0)*2D0 + 
     &             p1p2*(ampi**2 - p1p4*2D0) - 
     &             p1p4*2D0*(ame**2 - p2p3*1D0 + p2p4*2D0)) + 
     &          p1p4**2*p2p4*4D0) + 
     &       kp1*(ame**2*ampi**2*kp4 + ampi**2*kp4*p1p2 + 
     &          ame**2*kp4**2*2D0 + kp4**2*p1p2*2D0 - 
     &          ame**2*ampi**2*p2p4*2D0 - 
     &          ame**2*kp4*p2p4*2D0 - ampi**2*p1p2*p2p4*2D0 - 
     &          kp4*p1p2*p2p4*2D0 + kp4*p1p3*p2p4*2D0 - 
     &          p1p4*p2p3*p2p4*2D0 - p1p3*p2p4**2*2D0 + 
     &          (ame**2 + p1p2)*p2p4*p3p4*2D0 + 
     &          kp2**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &          kp3*((ame**2 + p1p2 + p1p4)*p2p4*2D0 - 
     &             (ame**2 + p1p2)*1D0*(ampi**2 + kp4*2D0)) + 
     &          kp2*(ampi**2*p2p3 - 
     &             p2p4*1D0*(ampi**2 + kp3*2D0 - kp4*2D0) + 
     &             p1p4*(ampi**2 + kp3*2D0 - kp4*2D0 + 
     &                p2p3*2D0) - 
     &             p1p3*1D0*(ampi**2 + p2p4*2D0)) - 
     &          kp4*p1p4*p2p4*4D0 + p1p4*p2p4**2*4D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D12ap0x12px35x34m1ppp*2D0*(m1**2 + p3p4*4D0)*
     &     (kp1**2*(p2p3**2 + (ampi**2 + kp3)*p2p4 + 
     &          kp2*(ampi**2 - p3p4*1D0) + 
     &          p2p3*(ampi**2 - kp3*1D0 - p2p4*1D0 - p3p4*2D0)
     &          ) + kp1*(kp3**2*(ame**2 + p1p2) + 
     &          ame**2*kp4*p2p3 + kp4*p1p2*p2p3 + 
     &          p2p3*p3p4*x12 - ame**2*ampi**2*kp4*1D0 - 
     &          ampi**2*kp4*p1p2*1D0 + 
     &          kp2**2*(p3p4 - ampi**2*1D0) + 
     &          kp2*(ampi**2*p2p4 - 
     &             p1p4*(ampi**2 + p2p3)*1D0 - 
     &             p2p3*1D0*(ampi**2 + kp4 - kp3*1D0) + 
     &             p1p3*(ampi**2 + kp4 + p2p4 - kp3*1D0)) - 
     &          ame**2*ampi**2*p2p3*2D0 - 
     &          ampi**2*p1p2*p2p3*2D0 + kp4*p1p3*p2p3*2D0 - 
     &          p1p4*p2p3**2*2D0 - p1p3*p2p3*p2p4*2D0 + 
     &          p1p3*p2p3**2*4D0 + 
     &          kp3*((ame**2 + p1p2)*(ampi**2 - kp4*1D0) - 
     &             p2p3*1D0*
     &              (ame**2 + p1p2 - p1p4*2D0 + p1p3*4D0))) - 
     &       kp2*1D0*(kp3**2*(ame**2 + p1p2) + 
     &          ampi**2*kp2*p1p3 + kp2*p1p3**2 + 
     &          ampi**2*kp2*p1p4 - kp2*p1p3*p1p4*1D0 - 
     &          ame**2*ampi**2*p1p3*2D0 - 
     &          ampi**2*p1p2*p1p3*2D0 - p1p3*p1p4*p2p3*2D0 - 
     &          p1p3**2*p2p4*2D0 + 
     &          p1p3*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 + 
     &          kp4*(-((ame**2 + p1p2)*1D0*
     &                (ampi**2 - p1p3*1D0)) + p1p3*p2p3*2D0)
     &           + p1p3**2*p2p3*4D0 + 
     &          kp3*(ampi**2*(ame**2 + p1p2) + kp2*p1p4 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             p1p3*1D0*
     &              (ame**2 + kp2 + p1p2 - p2p4*2D0 + 
     &                p2p3*4D0)))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C11ax12x340ppp*2D0*
     &     (kp1**2*(p2p3 - p2p4*1D0)*
     &        (p2p3 + p2p4 - kp3*1D0 - kp4*1D0 - 
     &          (ampi**2 + p3p4)*2D0) + 
     &       kp1*(kp3**2*(ame**2 + p1p2) - 
     &          kp4**2*(ame**2 + p1p2)*1D0 + 
     &          kp2*(p1p3*(ampi**2 + p2p4 + p3p4) - 
     &             p1p4*(ampi**2 + p2p3 + p3p4)*1D0 - 
     &             (ampi**2 + p3p4)*1D0*(p2p3 - p2p4*1D0))*2D0
     &            + (p2p3 + p2p4)*(p1p3*p2p3 - p1p4*p2p4*1D0)*
     &           4D0 + kp3*
     &           (kp2*(p2p3 + p2p4 - p1p3*1D0 - p1p4*1D0) + 
     &             ame**2*ampi**2*2D0 + ame**2*p3p4*2D0 - 
     &             p1p2*1D0*
     &              (p2p3 + p2p4 - (ampi**2 + p3p4)*2D0) - 
     &             (p2p3 + p2p4)*1D0*(ame**2 + p1p3*4D0)) + 
     &          kp4*(kp2*
     &              (p1p3 + p1p4 - p2p3*1D0 - p2p4*1D0) - 
     &             ame**2*ampi**2*2D0 - ame**2*p3p4*2D0 + 
     &             p1p2*
     &              (p2p3 + p2p4 - (ampi**2 + p3p4)*2D0) + 
     &             (p2p3 + p2p4)*(ame**2 + p1p4*4D0))) + 
     &       kp2*(kp4**2*(ame**2 + p1p2) + kp2*p1p4**2 - 
     &          kp3**2*(ame**2 + p1p2)*1D0 - 
     &          kp2*p1p3**2*1D0 + ampi**2*kp2*p1p3*2D0 - 
     &          ampi**2*kp2*p1p4*2D0 + 
     &          kp2*p3p4*(p1p3 - p1p4*1D0)*2D0 - 
     &          p1p3**2*p2p3*4D0 - p1p3*p1p4*p2p3*4D0 + 
     &          p1p3*p1p4*p2p4*4D0 + p1p4**2*p2p4*4D0 + 
     &          kp3*((ame**2 + kp2)*p1p3 + 
     &             p1p4*(ame**2 - kp2*1D0) - 
     &             ame**2*(ampi**2 + p3p4)*2D0 + 
     &             p1p2*
     &              (p1p3 + p1p4 - (ampi**2 + p3p4)*2D0) + 
     &             (p1p3 + p1p4)*p2p3*4D0) + 
     &          kp4*(ame**2*ampi**2*2D0 + ame**2*p3p4*2D0 - 
     &             p1p2*1D0*
     &              (p1p3 + p1p4 - (ampi**2 + p3p4)*2D0) + 
     &             p1p3*(kp2 - ame**2*1D0 - p2p4*4D0) - 
     &             p1p4*1D0*(ame**2 + kp2 + p2p4*4D0)))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D12ap0x12px45x34m1ppp*2D0*(m1**2 + p3p4*4D0)*
     &     (-(kp1**2*1D0*
     &          (p2p3*(ampi**2 + kp4 - p2p4*1D0) + 
     &            kp2*(ampi**2 - p3p4*1D0) + 
     &            p2p4*(ampi**2 + p2p4 - kp4*1D0 - p3p4*2D0)))
     &         + kp1*(ame**2*kp4*p2p4 + kp4*p1p2*p2p4 - 
     &          ame**2*ampi**2*kp4*1D0 - ame**2*kp4**2*1D0 - 
     &          ampi**2*kp4*p1p2*1D0 - kp4**2*p1p2*1D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0) + 
     &          kp2*(p1p3*(ampi**2 + p2p4) - 
     &             ampi**2*p2p3*1D0 + 
     &             p2p4*(ampi**2 + kp3 - kp4*1D0) - 
     &             p1p4*1D0*(ampi**2 + kp3 + p2p3 - kp4*1D0))
     &           + ame**2*ampi**2*p2p4*2D0 + 
     &          ampi**2*p1p2*p2p4*2D0 - kp4*p1p3*p2p4*2D0 + 
     &          p1p4*p2p3*p2p4*2D0 + p1p3*p2p4**2*2D0 - 
     &          (ame**2 + p1p2)*p2p4*p3p4*2D0 + 
     &          kp3*((ampi**2 + kp4)*(ame**2 + p1p2) - 
     &             p2p4*1D0*(ame**2 + p1p2 + p1p4*2D0)) + 
     &          kp4*p1p4*p2p4*4D0 - p1p4*p2p4**2*4D0) + 
     &       kp2*(kp4**2*(ame**2 + p1p2) + ampi**2*kp2*p1p3 + 
     &          ampi**2*kp2*p1p4 + kp2*p1p4**2 - 
     &          kp2*p1p3*p1p4*1D0 - ame**2*ampi**2*p1p4*2D0 - 
     &          ampi**2*p1p2*p1p4*2D0 - p1p4**2*p2p3*2D0 - 
     &          p1p3*p1p4*p2p4*2D0 + 
     &          p1p4*p3p4*(ame**2 + p1p2 - kp2*1D0)*2D0 + 
     &          kp3*(-((ame**2 + p1p2)*1D0*
     &                (ampi**2 + kp4 - p1p4*1D0)) + 
     &             p1p4*p2p4*2D0) + p1p4**2*p2p4*4D0 + 
     &          kp4*(ame**2*ampi**2 + 
     &             p1p2*(ampi**2 - p1p4*1D0) + 
     &             kp2*(p1p3 - p1p4*1D0) - 
     &             p1p4*1D0*(ame**2 - p2p3*2D0 + p2p4*4D0)))))
     &    /(kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C2apx35x12pm1p*(-(kp1**2*1D0*
     &          (p2p3**2*3D0 + kp2*(ampi**2 - p3p4*1D0)*3D0 + 
     &            p2p4*(ampi**2 + kp3*3D0) + 
     &            p2p3*(ampi**2*2D0 - kp3*3D0 - p2p4*3D0 - 
     &               p3p4*3D0))) + 
     &       kp2*(ampi**2*kp2*p1p4 + ampi**2*kp2*p1p3*2D0 + 
     &          p1p3*p1p4*p2p4*2D0 + 
     &          kp3**2*(ame**2 + p1p2)*3D0 + 
     &          kp2*p1p3**2*3D0 - kp2*p1p3*p1p4*3D0 - 
     &          kp4*1D0*(ame**2*ampi**2 + 
     &             p1p2*(ampi**2 - p1p3*3D0) + 
     &             p1p3*(p2p4 - (ame**2 + p2p3)*3D0)) - 
     &          ame**2*ampi**2*p1p3*4D0 - 
     &          ampi**2*p1p2*p1p3*4D0 - p1p3*p1p4*p2p3*4D0 - 
     &          p1p3**2*p2p4*4D0 + 
     &          p1p3*p3p4*
     &           (-(kp2*3D0) + (ame**2 + p1p2)*4D0) + 
     &          kp3*(-((ame**2 + p1p2)*p3p4*1D0) + 
     &             ampi**2*(ame**2 + p1p2)*2D0 - 
     &             kp4*(ame**2 + p1p2)*3D0 + 
     &             p1p4*(p2p3 - p2p4*1D0 + kp2*3D0) + 
     &             p1p3*(-((ame**2 + kp2 + p1p2 + p2p3*2D0)*
     &                   3D0) + p2p4*4D0)) + p1p3**2*p2p3*6D0)
     &         + kp1*(ame**2*ampi**2*kp4 + ampi**2*kp4*p1p2 + 
     &          kp4*p1p4*p2p3 - p1p4*p2p3*p2p4*2D0 - 
     &          kp3**2*(ame**2 + p1p2)*3D0 - 
     &          ame**2*kp4*p2p3*3D0 - kp4*p1p2*p2p3*3D0 - 
     &          kp4*p1p3*p2p3*3D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0)*3D0 + 
     &          kp2*(ampi**2*p2p3 + p2p4*p3p4 - 
     &             ampi**2*p2p4*2D0 - kp3*p2p3*3D0 + 
     &             kp4*p2p3*3D0 + 
     &             p1p4*(-(p3p4*1D0) + ampi**2*2D0 + 
     &                p2p3*3D0) - 
     &             p1p3*1D0*
     &              (ampi**2 - kp3*3D0 + kp4*3D0 + p2p4*3D0))
     &           + ame**2*ampi**2*p2p3*4D0 + 
     &          ampi**2*p1p2*p2p3*4D0 + p1p4*p2p3**2*4D0 + 
     &          p1p3*p2p3*p2p4*4D0 - 
     &          (ame**2 + p1p2)*p2p3*p3p4*4D0 + 
     &          kp3*((ame**2 + p1p2)*p3p4 + 
     &             p2p4*(p1p4 - p1p3*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             kp4*(ame**2 + p1p2)*3D0 + 
     &             p2p3*((ame**2 + p1p2 + p1p3*2D0)*3D0 - 
     &                p1p4*4D0)) - p1p3*p2p3**2*6D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D2ap0x12px35x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (kp1**2*(p2p3**2*3D0 + 
     &          kp2*(ampi**2 - p3p4*1D0)*3D0 + 
     &          p2p4*(ampi**2 + kp3*3D0) + 
     &          p2p3*(ampi**2*2D0 - kp3*3D0 - p2p4*3D0 - 
     &             p3p4*3D0)) + 
     &       kp2*(-(ampi**2*kp2*p1p4*1D0) - 
     &          ampi**2*kp2*p1p3*2D0 - p1p3*p1p4*p2p4*2D0 - 
     &          kp3**2*(ame**2 + p1p2)*3D0 - 
     &          kp2*p1p3**2*3D0 + kp2*p1p3*p1p4*3D0 + 
     &          kp4*(ame**2*ampi**2 + 
     &             p1p2*(ampi**2 - p1p3*3D0) + 
     &             p1p3*(p2p4 - (ame**2 + p2p3)*3D0)) + 
     &          ame**2*ampi**2*p1p3*4D0 + 
     &          ampi**2*p1p2*p1p3*4D0 + p1p3*p1p4*p2p3*4D0 + 
     &          p1p3**2*p2p4*4D0 + 
     &          p1p3*p3p4*(kp2*3D0 - (ame**2 + p1p2)*4D0) + 
     &          kp3*((ame**2 + p1p2)*p3p4 - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             kp4*(ame**2 + p1p2)*3D0 + 
     &             p1p4*(p2p4 - p2p3*1D0 - kp2*3D0) + 
     &             p1p3*((ame**2 + kp2 + p1p2 + p2p3*2D0)*
     &                 3D0 - p2p4*4D0)) - p1p3**2*p2p3*6D0) - 
     &       kp1*1D0*(ame**2*ampi**2*kp4 + ampi**2*kp4*p1p2 + 
     &          kp4*p1p4*p2p3 - p1p4*p2p3*p2p4*2D0 - 
     &          kp3**2*(ame**2 + p1p2)*3D0 - 
     &          ame**2*kp4*p2p3*3D0 - kp4*p1p2*p2p3*3D0 - 
     &          kp4*p1p3*p2p3*3D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0)*3D0 + 
     &          kp2*(ampi**2*p2p3 + p2p4*p3p4 - 
     &             ampi**2*p2p4*2D0 - kp3*p2p3*3D0 + 
     &             kp4*p2p3*3D0 + 
     &             p1p4*(-(p3p4*1D0) + ampi**2*2D0 + 
     &                p2p3*3D0) - 
     &             p1p3*1D0*
     &              (ampi**2 - kp3*3D0 + kp4*3D0 + p2p4*3D0))
     &           + ame**2*ampi**2*p2p3*4D0 + 
     &          ampi**2*p1p2*p2p3*4D0 + p1p4*p2p3**2*4D0 + 
     &          p1p3*p2p3*p2p4*4D0 - 
     &          (ame**2 + p1p2)*p2p3*p3p4*4D0 + 
     &          kp3*((ame**2 + p1p2)*p3p4 + 
     &             p2p4*(p1p4 - p1p3*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             kp4*(ame**2 + p1p2)*3D0 + 
     &             p2p3*((ame**2 + p1p2 + p1p3*2D0)*3D0 - 
     &                p1p4*4D0)) - p1p3*p2p3**2*6D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D3ap0x12px45x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (-(kp1**2*1D0*
     &          (p2p3**2 + (ampi**2 + kp3)*p2p4 + 
     &            kp2*(ampi**2 - p3p4*1D0) - 
     &            p2p3*1D0*
     &             (kp3 + p2p4 - p3p4*1D0 + ampi**2*2D0))) + 
     &       kp2*(kp3**2*(ame**2 + p1p2) + kp2*p1p3**2 + 
     &          ampi**2*kp2*p1p4 + kp2*p1p3*p3p4 + 
     &          p1p4*p3p4*x12 - kp2*p1p3*p1p4*1D0 + 
     &          ame**2*ampi**2*p1p3*2D0 - 
     &          ampi**2*kp2*p1p3*2D0 + 
     &          ampi**2*p1p2*p1p3*2D0 - 
     &          ame**2*ampi**2*p1p4*2D0 - 
     &          ampi**2*p1p2*p1p4*2D0 - p1p3**2*p2p3*2D0 - 
     &          p1p4**2*p2p3*2D0 + p1p3**2*p2p4*2D0 - 
     &          ame**2*p1p3*p3p4*2D0 - p1p2*p1p3*p3p4*2D0 + 
     &          kp3*((ame**2 + p1p2)*p3p4 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             p1p3*1D0*
     &              (ame**2 + kp2 + p1p2 - p2p3*2D0 + 
     &                p2p4*2D0) + p1p4*(kp2 + p2p4 - p2p3*3D0)
     &             ) + kp4*
     &           (p1p4*p2p3*2D0 + 
     &             ame**2*(ampi**2 - p3p4*2D0) + 
     &             p1p2*(ampi**2 + p1p3 - p3p4*2D0) + 
     &             p1p3*(ame**2 - p2p3*3D0 + p2p4*3D0)) - 
     &          p1p3*p1p4*p2p4*4D0 + p1p3*p1p4*p2p3*6D0) + 
     &       kp1*(-(ame**2*ampi**2*kp4*1D0) - 
     &          ampi**2*kp4*p1p2*1D0 - 
     &          kp3**2*(ame**2 + p1p2)*1D0 - 
     &          ame**2*kp4*p2p3*1D0 - kp4*p1p2*p2p3*1D0 + 
     &          p3p4*x12*(kp4 + p2p3 - p2p4*1D0) + 
     &          kp2**2*(ampi**2 - p3p4*1D0) - 
     &          ame**2*ampi**2*p2p3*2D0 - 
     &          ampi**2*p1p2*p2p3*2D0 + p1p3*p2p3**2*2D0 - 
     &          p1p4*p2p3**2*2D0 + ame**2*ampi**2*p2p4*2D0 + 
     &          ampi**2*p1p2*p2p4*2D0 - kp4*p1p3*p2p4*2D0 + 
     &          p1p3*p2p4**2*2D0 + 
     &          kp2*(p3p4*(p1p4 - p2p4*1D0) + 
     &             p1p3*(ampi**2 + kp3 - kp4*1D0 - p2p4*1D0 - 
     &                p3p4*2D0) + 
     &             p2p3*(kp4 + p1p4 - ampi**2*1D0 - kp3*1D0 + 
     &                p3p4*2D0)) + kp4*p1p3*p2p3*3D0 - 
     &          kp4*p1p4*p2p3*3D0 + 
     &          kp3*(kp4*(ame**2 + p1p2) - 
     &             (ame**2 + p1p2)*p3p4*1D0 + 
     &             p2p3*(ame**2 + p1p2 - p1p3*2D0 + 
     &                p1p4*2D0) + 
     &             p2p4*(-(p1p4*1D0) + p1p3*3D0)) + 
     &          p1p4*p2p3*p2p4*4D0 - p1p3*p2p3*p2p4*6D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C2ax45px12pm1p*(kp1**2*
     &        (p2p3**2 + (ampi**2 + kp3)*p2p4 + 
     &          kp2*(ampi**2 - p3p4*1D0) - 
     &          p2p3*1D0*(kp3 + p2p4 - p3p4*1D0 + ampi**2*2D0)
     &          ) + kp2*(kp2*p1p3*p1p4 - 
     &          kp3**2*(ame**2 + p1p2)*1D0 - 
     &          kp2*p1p3**2*1D0 - ampi**2*kp2*p1p4*1D0 - 
     &          ame**2*ampi**2*p1p3*2D0 + 
     &          ampi**2*kp2*p1p3*2D0 - 
     &          ampi**2*p1p2*p1p3*2D0 + 
     &          ame**2*ampi**2*p1p4*2D0 + 
     &          ampi**2*p1p2*p1p4*2D0 + p1p3**2*p2p3*2D0 + 
     &          p1p4**2*p2p3*2D0 - p1p3**2*p2p4*2D0 + 
     &          p3p4*(-((ame**2 + p1p2)*p1p4*2D0) + 
     &             p1p3*(-(kp2*1D0) + (ame**2 + p1p2)*2D0)) + 
     &          kp4*(-(p1p4*p2p3*2D0) - 
     &             ame**2*1D0*(ampi**2 + p1p3 - p3p4*2D0) - 
     &             p1p2*1D0*(ampi**2 + p1p3 - p3p4*2D0) + 
     &             p1p3*(p2p3 - p2p4*1D0)*3D0) + 
     &          kp3*(kp4*(ame**2 + p1p2) - 
     &             (ame**2 + p1p2)*p3p4*1D0 + 
     &             p1p3*(ame**2 + kp2 + p1p2 - p2p3*2D0 + 
     &                p2p4*2D0) - 
     &             p1p4*1D0*(kp2 + p2p4 - p2p3*3D0)) + 
     &          p1p3*p1p4*p2p4*4D0 - p1p3*p1p4*p2p3*6D0) + 
     &       kp1*(ame**2*ampi**2*kp4 + ampi**2*kp4*p1p2 + 
     &          kp3**2*(ame**2 + p1p2) + ame**2*kp4*p2p3 + 
     &          kp4*p1p2*p2p3 + kp2**2*(p3p4 - ampi**2*1D0) + 
     &          ame**2*ampi**2*p2p3*2D0 + 
     &          ampi**2*p1p2*p2p3*2D0 - p1p3*p2p3**2*2D0 + 
     &          p1p4*p2p3**2*2D0 - ame**2*ampi**2*p2p4*2D0 - 
     &          ampi**2*p1p2*p2p4*2D0 + kp4*p1p3*p2p4*2D0 - 
     &          p1p3*p2p4**2*2D0 - 
     &          (ame**2 + p1p2)*p3p4*(kp4 + p2p3 - p2p4*1D0)*
     &           2D0 + kp2*
     &           (p3p4*(p2p4 - p1p4*1D0) + 
     &             p2p3*(ampi**2 + kp3 - kp4*1D0 - p1p4*1D0 - 
     &                p3p4*2D0) + 
     &             p1p3*(kp4 + p2p4 - ampi**2*1D0 - kp3*1D0 + 
     &                p3p4*2D0)) - kp4*p1p3*p2p3*3D0 + 
     &          kp4*p1p4*p2p3*3D0 + 
     &          kp3*((ame**2 + p1p2)*p3p4 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             p2p3*1D0*
     &              (ame**2 + p1p2 - p1p3*2D0 + p1p4*2D0) + 
     &             p2p4*(p1p4 - p1p3*3D0)) - 
     &          p1p4*p2p3*p2p4*4D0 + p1p3*p2p3*p2p4*6D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C2ax35px12pm1p*(-(kp1**2*1D0*
     &          (p2p3*(ampi**2 + kp4 - p2p4*1D0) + 
     &            kp2*(ampi**2 - p3p4*1D0) + 
     &            p2p4*(p2p4 + p3p4 - kp4*1D0 - ampi**2*2D0)))
     &         + kp2*(kp4**2*(ame**2 + p1p2) + 
     &          ampi**2*kp2*p1p3 + kp2*p1p4**2 - 
     &          kp2*p1p3*p1p4*1D0 - ame**2*ampi**2*p1p3*2D0 - 
     &          ampi**2*p1p2*p1p3*2D0 + 
     &          ame**2*ampi**2*p1p4*2D0 - 
     &          ampi**2*kp2*p1p4*2D0 + 
     &          ampi**2*p1p2*p1p4*2D0 + p1p4**2*p2p3*2D0 - 
     &          p1p3**2*p2p4*2D0 - p1p4**2*p2p4*2D0 + 
     &          p3p4*((ame**2 + p1p2)*p1p3*2D0 + 
     &             p1p4*(kp2 - (ame**2 + p1p2)*2D0)) + 
     &          kp4*((ame**2 + p1p2)*p3p4 + 
     &             kp2*(p1p3 - p1p4*1D0) - 
     &             p1p4*1D0*
     &              (ame**2 + p1p2 + p2p3*2D0 - p2p4*2D0) + 
     &             p1p3*(p2p3 - p2p4*3D0)) + 
     &          kp3*(ampi**2*(ame**2 + p1p2) - 
     &             kp4*(ame**2 + p1p2)*1D0 + p1p3*p2p4*2D0 - 
     &             (ame**2 + p1p2)*p3p4*2D0 + 
     &             p1p4*(ame**2 + p1p2 + p2p3*3D0 - p2p4*3D0))
     &            - p1p3*p1p4*p2p3*4D0 + p1p3*p1p4*p2p4*6D0)
     &        + kp1*(ame**2*kp4*p2p4 + kp4*p1p2*p2p4 - 
     &          ame**2*kp4**2*1D0 - kp4**2*p1p2*1D0 - 
     &          kp4*p1p3*p2p3*1D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0) + 
     &          ame**2*ampi**2*p2p3*2D0 + 
     &          ampi**2*p1p2*p2p3*2D0 + p1p4*p2p3**2*2D0 - 
     &          ame**2*ampi**2*p2p4*2D0 - 
     &          ampi**2*p1p2*p2p4*2D0 + kp4*p1p3*p2p4*2D0 - 
     &          kp4*p1p4*p2p4*2D0 - p1p3*p2p4**2*2D0 + 
     &          p1p4*p2p4**2*2D0 - 
     &          (ame**2 + p1p2)*p3p4*1D0*
     &           (kp4 + p2p3*2D0 - p2p4*2D0) + 
     &          kp2*(p2p4*
     &              (kp3 + p1p3 - ampi**2*1D0 - kp4*1D0) + 
     &             p3p4*(p1p3 - p2p3*1D0 + p2p4*2D0) + 
     &             p1p4*(ampi**2 + kp4 - kp3*1D0 - p2p3*1D0 - 
     &                p3p4*2D0)) + kp4*p1p4*p2p3*3D0 + 
     &          kp3*(kp4*(ame**2 + p1p2) - 
     &             ampi**2*(ame**2 + p1p2)*1D0 - 
     &             p1p4*p2p3*2D0 + (ame**2 + p1p2)*p3p4*2D0 - 
     &             p2p4*1D0*
     &              (ame**2 + p1p2 + p1p3*3D0 - p1p4*3D0)) + 
     &          p1p3*p2p3*p2p4*4D0 - p1p4*p2p3*p2p4*6D0)))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D3ap0x12px35x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (kp1**2*(p2p3*(ampi**2 + kp4 - p2p4*1D0) + 
     &          kp2*(ampi**2 - p3p4*1D0) + 
     &          p2p4*(p2p4 + p3p4 - kp4*1D0 - ampi**2*2D0)) - 
     &       kp2*1D0*(kp4**2*(ame**2 + p1p2) + 
     &          ampi**2*kp2*p1p3 + kp2*p1p4**2 + 
     &          kp2*p1p4*p3p4 + p1p3*p3p4*x12 - 
     &          kp2*p1p3*p1p4*1D0 - ame**2*ampi**2*p1p3*2D0 - 
     &          ampi**2*p1p2*p1p3*2D0 + 
     &          ame**2*ampi**2*p1p4*2D0 - 
     &          ampi**2*kp2*p1p4*2D0 + 
     &          ampi**2*p1p2*p1p4*2D0 + p1p4**2*p2p3*2D0 - 
     &          p1p3**2*p2p4*2D0 - p1p4**2*p2p4*2D0 - 
     &          ame**2*p1p4*p3p4*2D0 - p1p2*p1p4*p3p4*2D0 + 
     &          kp4*((ame**2 + p1p2)*p3p4 + 
     &             kp2*(p1p3 - p1p4*1D0) - 
     &             p1p4*1D0*
     &              (ame**2 + p1p2 + p2p3*2D0 - p2p4*2D0) + 
     &             p1p3*(p2p3 - p2p4*3D0)) + 
     &          kp3*(ampi**2*(ame**2 + p1p2) - 
     &             kp4*(ame**2 + p1p2)*1D0 + p1p3*p2p4*2D0 - 
     &             (ame**2 + p1p2)*p3p4*2D0 + 
     &             p1p4*(ame**2 + p1p2 + p2p3*3D0 - p2p4*3D0))
     &            - p1p3*p1p4*p2p3*4D0 + p1p3*p1p4*p2p4*6D0)
     &        + kp1*(ame**2*kp4**2 + kp4**2*p1p2 + 
     &          kp4*p1p3*p2p3 - ame**2*kp4*p2p4*1D0 - 
     &          kp4*p1p2*p2p4*1D0 + 
     &          kp2**2*(p3p4 - ampi**2*1D0) - 
     &          ame**2*ampi**2*p2p3*2D0 - 
     &          ampi**2*p1p2*p2p3*2D0 - p1p4*p2p3**2*2D0 + 
     &          ame**2*ampi**2*p2p4*2D0 + 
     &          ampi**2*p1p2*p2p4*2D0 - kp4*p1p3*p2p4*2D0 + 
     &          kp4*p1p4*p2p4*2D0 + p1p3*p2p4**2*2D0 - 
     &          p1p4*p2p4**2*2D0 + 
     &          (ame**2 + p1p2)*p3p4*
     &           (kp4 + p2p3*2D0 - p2p4*2D0) + 
     &          kp2*(p2p4*
     &              (ampi**2 + kp4 - kp3*1D0 - p1p3*1D0) + 
     &             p3p4*(p2p3 - p1p3*1D0 - p2p4*2D0) + 
     &             p1p4*(kp3 + p2p3 - ampi**2*1D0 - kp4*1D0 + 
     &                p3p4*2D0)) - kp4*p1p4*p2p3*3D0 + 
     &          kp3*(ampi**2*(ame**2 + p1p2) - 
     &             kp4*(ame**2 + p1p2)*1D0 + p1p4*p2p3*2D0 - 
     &             (ame**2 + p1p2)*p3p4*2D0 + 
     &             p2p4*(ame**2 + p1p2 + p1p3*3D0 - p1p4*3D0))
     &            - p1p3*p2p3*p2p4*4D0 + p1p4*p2p3*p2p4*6D0)))
     &    /(kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D1ap0x12px35x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (kp1**2*(p2p3**2 + (ampi**2 + kp3)*p2p4 + 
     &          kp2*(ampi**2 - p3p4*1D0) + 
     &          p2p3*(-(kp3*1D0) - p2p4*1D0 + ampi**2*2D0 - 
     &             p3p4*3D0)) - 
     &       kp1*1D0*(ame**2*ampi**2*kp4 + ampi**2*kp4*p1p2 + 
     &          kp4*p1p4*p2p3 - kp3**2*(ame**2 + p1p2)*1D0 - 
     &          ame**2*kp4*p2p3*1D0 - kp4*p1p2*p2p3*1D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0) - 
     &          p1p4*p2p3*p2p4*2D0 + 
     &          kp2*(ampi**2*p2p3 + kp4*p2p3 + p2p4*p3p4 - 
     &             kp3*p2p3*1D0 - 
     &             p1p3*1D0*
     &              (ampi**2 + kp4 + p2p4 - kp3*1D0) - 
     &             ampi**2*p2p4*2D0 + 
     &             p1p4*(p2p3 - p3p4*1D0 + ampi**2*2D0)) - 
     &          kp4*p1p3*p2p3*3D0 + ame**2*ampi**2*p2p3*4D0 + 
     &          ampi**2*p1p2*p2p3*4D0 + p1p4*p2p3**2*4D0 + 
     &          p1p3*p2p3*p2p4*4D0 - 
     &          (ame**2 + p1p2)*p2p3*p3p4*4D0 - 
     &          p1p3*p2p3**2*6D0 + 
     &          kp3*(kp4*(ame**2 + p1p2) + 
     &             (ame**2 + p1p2)*p3p4 + 
     &             p2p4*(p1p4 - p1p3*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             p2p3*(ame**2 + p1p2 - p1p4*4D0 + p1p3*6D0))
     &          ) + kp2*(kp2*p1p3*p1p4 - 
     &          kp3**2*(ame**2 + p1p2)*1D0 - 
     &          kp2*p1p3**2*1D0 - ampi**2*kp2*p1p4*1D0 - 
     &          ampi**2*kp2*p1p3*2D0 - p1p3*p1p4*p2p4*2D0 + 
     &          kp4*(ame**2*ampi**2 + 
     &             p1p2*(ampi**2 - p1p3*1D0) + 
     &             p1p3*(p2p4 - ame**2*1D0 - p2p3*3D0)) + 
     &          ame**2*ampi**2*p1p3*4D0 + 
     &          ampi**2*p1p2*p1p3*4D0 + p1p3*p1p4*p2p3*4D0 + 
     &          p1p3**2*p2p4*4D0 + 
     &          p1p3*p3p4*(kp2*3D0 - (ame**2 + p1p2)*4D0) - 
     &          p1p3**2*p2p3*6D0 + 
     &          kp3*(kp4*(ame**2 + p1p2) + 
     &             (ame**2 + p1p2)*p3p4 - 
     &             p1p4*1D0*(kp2 + p2p3 - p2p4*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             p1p3*(ame**2 + kp2 + p1p2 - p2p4*4D0 + 
     &                p2p3*6D0)))))/
     &   (kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D2ap0x12px45x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (-(kp1**2*1D0*
     &          (kp2*(ampi**2 - p3p4*1D0)*3D0 + 
     &            p2p3*(ampi**2 + kp4*3D0 - p2p4*3D0) + 
     &            p2p4*(ampi**2*2D0 - kp4*3D0 + p2p4*3D0 - 
     &               p3p4*3D0))) + 
     &       kp1*(kp4*p1p3*p2p3 - kp4*p1p4*p2p3*1D0 - 
     &          ame**2*ampi**2*kp4*2D0 - 
     &          ampi**2*kp4*p1p2*2D0 - p1p3*p2p3*p2p4*2D0 - 
     &          ame**2*kp4**2*3D0 - kp4**2*p1p2*3D0 + 
     &          ame**2*kp4*p2p4*3D0 + kp4*p1p2*p2p4*3D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0)*3D0 + 
     &          kp3*((ame**2 + p1p2)*(ampi**2 + kp4*3D0) + 
     &             p2p4*(p1p3 - p1p2*3D0 - 
     &                (ame**2 + p1p4)*3D0)) + 
     &          kp2*(ampi**2*p2p4 + p2p3*p3p4 - 
     &             ampi**2*p2p3*2D0 + kp3*p2p4*3D0 - 
     &             kp4*p2p4*3D0 - 
     &             p1p4*1D0*
     &              (ampi**2 + kp3*3D0 - kp4*3D0 + p2p3*3D0)
     &              + p1p3*
     &              (-(p3p4*1D0) + ampi**2*2D0 + p2p4*3D0)) + 
     &          ame**2*ampi**2*p2p4*4D0 + 
     &          ampi**2*p1p2*p2p4*4D0 - kp4*p1p3*p2p4*4D0 + 
     &          p1p4*p2p3*p2p4*4D0 + p1p3*p2p4**2*4D0 + 
     &          (ame**2 + p1p2)*p3p4*(kp4 - p2p4*4D0) + 
     &          kp4*p1p4*p2p4*6D0 - p1p4*p2p4**2*6D0) + 
     &       kp2*(ampi**2*kp2*p1p3 + ampi**2*kp2*p1p4*2D0 + 
     &          p1p3*p1p4*p2p3*2D0 + 
     &          kp4**2*(ame**2 + p1p2)*3D0 - 
     &          kp2*p1p3*p1p4*3D0 + kp2*p1p4**2*3D0 - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) + 
     &             kp4*(ame**2 + p1p2)*3D0 + 
     &             p1p4*(p2p3 - p1p2*3D0 - 
     &                (ame**2 + p2p4)*3D0)) - 
     &          ame**2*ampi**2*p1p4*4D0 - 
     &          ampi**2*p1p2*p1p4*4D0 - p1p4**2*p2p3*4D0 - 
     &          p1p3*p1p4*p2p4*4D0 + 
     &          p1p4*p3p4*
     &           (-(kp2*3D0) + (ame**2 + p1p2)*4D0) + 
     &          p1p4**2*p2p4*6D0 + 
     &          kp4*(-(ame**2*p3p4*1D0) + 
     &             p1p3*(p2p4 - p2p3*1D0) + 
     &             ame**2*ampi**2*2D0 + 
     &             kp2*(p1p3 - p1p4*1D0)*3D0 + 
     &             p1p2*(-(p3p4*1D0) + ampi**2*2D0 - 
     &                p1p4*3D0) + 
     &             p1p4*(-(ame**2*3D0) + p2p3*4D0 - p2p4*6D0))
     &          )))/(kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (D1ap0x12px45x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (-(kp1**2*1D0*
     &          (p2p3*(ampi**2 + kp4 - p2p4*1D0) + 
     &            kp2*(ampi**2 - p3p4*1D0) + 
     &            p2p4*(p2p4 - kp4*1D0 + ampi**2*2D0 - 
     &               p3p4*3D0))) + 
     &       kp1*(kp4*p1p3*p2p3 + ame**2*kp4*p2p4 + 
     &          kp4*p1p2*p2p4 - ame**2*kp4**2*1D0 - 
     &          kp4**2*p1p2*1D0 - kp4*p1p4*p2p3*1D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0) - 
     &          ame**2*ampi**2*kp4*2D0 - 
     &          ampi**2*kp4*p1p2*2D0 - p1p3*p2p3*p2p4*2D0 + 
     &          kp2*(ampi**2*p2p4 + kp3*p2p4 + p2p3*p3p4 - 
     &             kp4*p2p4*1D0 - 
     &             p1p4*1D0*
     &              (ampi**2 + kp3 + p2p3 - kp4*1D0) - 
     &             ampi**2*p2p3*2D0 + 
     &             p1p3*(p2p4 - p3p4*1D0 + ampi**2*2D0)) + 
     &          kp3*((ampi**2 + kp4)*(ame**2 + p1p2) - 
     &             p2p4*1D0*
     &              (ame**2 + p1p2 - p1p3*1D0 + p1p4*3D0)) + 
     &          ame**2*ampi**2*p2p4*4D0 + 
     &          ampi**2*p1p2*p2p4*4D0 - kp4*p1p3*p2p4*4D0 + 
     &          p1p4*p2p3*p2p4*4D0 + p1p3*p2p4**2*4D0 + 
     &          (ame**2 + p1p2)*p3p4*(kp4 - p2p4*4D0) + 
     &          kp4*p1p4*p2p4*6D0 - p1p4*p2p4**2*6D0) + 
     &       kp2*(kp4**2*(ame**2 + p1p2) + ampi**2*kp2*p1p3 + 
     &          kp2*p1p4**2 - kp2*p1p3*p1p4*1D0 + 
     &          ampi**2*kp2*p1p4*2D0 + p1p3*p1p4*p2p3*2D0 - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) + 
     &             kp4*(ame**2 + p1p2) - 
     &             p1p4*1D0*
     &              (ame**2 + p1p2 - p2p3*1D0 + p2p4*3D0)) - 
     &          ame**2*ampi**2*p1p4*4D0 - 
     &          ampi**2*p1p2*p1p4*4D0 - p1p4**2*p2p3*4D0 - 
     &          p1p3*p1p4*p2p4*4D0 + 
     &          p1p4*p3p4*
     &           (-(kp2*3D0) + (ame**2 + p1p2)*4D0) + 
     &          p1p4**2*p2p4*6D0 + 
     &          kp4*(-(ame**2*p3p4*1D0) + 
     &             kp2*(p1p3 - p1p4*1D0) + 
     &             p1p3*(p2p4 - p2p3*1D0) + 
     &             ame**2*ampi**2*2D0 + 
     &             p1p2*(-(p1p4*1D0) - p3p4*1D0 + 
     &                ampi**2*2D0) - 
     &             p1p4*1D0*(ame**2 - p2p3*4D0 + p2p4*6D0)))))
     &    /(kp1*kp2*(ampi**2 + p3p4)*x12) + 
     &  (C2apx45x12pm1p*(kp1**2*
     &        (kp2*(ampi**2 - p3p4*1D0)*3D0 + 
     &          p2p3*(ampi**2 + kp4*3D0 - p2p4*3D0) + 
     &          p2p4*(ampi**2*2D0 - kp4*3D0 + p2p4*3D0 - 
     &             p3p4*3D0)) - 
     &       kp1*1D0*(kp4*p1p3*p2p3 - kp4*p1p4*p2p3*1D0 - 
     &          ame**2*ampi**2*kp4*2D0 - 
     &          ampi**2*kp4*p1p2*2D0 - p1p3*p2p3*p2p4*2D0 - 
     &          ame**2*kp4**2*3D0 - kp4**2*p1p2*3D0 + 
     &          ame**2*kp4*p2p4*3D0 + kp4*p1p2*p2p4*3D0 + 
     &          kp2**2*(ampi**2 - p3p4*1D0)*3D0 + 
     &          kp3*((ame**2 + p1p2)*(ampi**2 + kp4*3D0) + 
     &             p2p4*(p1p3 - p1p2*3D0 - 
     &                (ame**2 + p1p4)*3D0)) + 
     &          kp2*(ampi**2*p2p4 + p2p3*p3p4 - 
     &             ampi**2*p2p3*2D0 + kp3*p2p4*3D0 - 
     &             kp4*p2p4*3D0 - 
     &             p1p4*1D0*
     &              (ampi**2 + kp3*3D0 - kp4*3D0 + p2p3*3D0)
     &              + p1p3*
     &              (-(p3p4*1D0) + ampi**2*2D0 + p2p4*3D0)) + 
     &          ame**2*ampi**2*p2p4*4D0 + 
     &          ampi**2*p1p2*p2p4*4D0 - kp4*p1p3*p2p4*4D0 + 
     &          p1p4*p2p3*p2p4*4D0 + p1p3*p2p4**2*4D0 + 
     &          (ame**2 + p1p2)*p3p4*(kp4 - p2p4*4D0) + 
     &          kp4*p1p4*p2p4*6D0 - p1p4*p2p4**2*6D0) + 
     &       kp2*(-(ampi**2*kp2*p1p3*1D0) - 
     &          ampi**2*kp2*p1p4*2D0 - p1p3*p1p4*p2p3*2D0 - 
     &          kp4**2*(ame**2 + p1p2)*3D0 + 
     &          kp2*p1p3*p1p4*3D0 - kp2*p1p4**2*3D0 + 
     &          kp3*(ampi**2*(ame**2 + p1p2) + 
     &             kp4*(ame**2 + p1p2)*3D0 + 
     &             p1p4*(p2p3 - p1p2*3D0 - 
     &                (ame**2 + p2p4)*3D0)) + 
     &          ame**2*ampi**2*p1p4*4D0 + 
     &          ampi**2*p1p2*p1p4*4D0 + p1p4**2*p2p3*4D0 + 
     &          p1p3*p1p4*p2p4*4D0 + 
     &          p1p4*p3p4*(kp2*3D0 - (ame**2 + p1p2)*4D0) - 
     &          p1p4**2*p2p4*6D0 + 
     &          kp4*(ame**2*p3p4 + p1p3*(p2p3 - p2p4*1D0) - 
     &             ame**2*ampi**2*2D0 + 
     &             kp2*(p1p4 - p1p3*1D0)*3D0 + 
     &             p1p2*(p3p4 - ampi**2*2D0 + p1p4*3D0) + 
     &             p1p4*(ame**2*3D0 - p2p3*4D0 + p2p4*6D0)))))
     &    /(kp1*kp2*(ampi**2 + p3p4)*x12)
