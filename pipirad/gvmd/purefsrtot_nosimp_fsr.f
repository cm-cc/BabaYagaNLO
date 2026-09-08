        elmat2 = -((C0apx35x12pm1p*1D0*
     &       (ampi**2*kp1*kp2*kp4 + ame**2*kp3**2*kp4 + 
     &         ame**2*ampi**2*kp4**2 + kp3**2*kp4*p1p2 + 
     &         ampi**2*kp4**2*p1p2 + ampi**2*kp2*kp4*p1p3 + 
     &         ampi**2*kp1*kp4*p2p3 + ampi**2*kp4*p1p3*p2p3 + 
     &         kp3*kp4*p1p4*p2p3 + kp3*kp4*p1p3*p2p4 + 
     &         ampi**2*kp4*p1p4*p2p4 + 
     &         ame**2*ampi**2*kp3*p3p4 + kp1*kp2*kp3*p3p4 + 
     &         ame**2*ampi**2*kp4*p3p4 + 
     &         ampi**2*kp3*p1p2*p3p4 + 
     &         ampi**2*kp4*p1p2*p3p4 + kp3*p1p4*p2p3*p3p4 + 
     &         kp3*p1p3*p2p4*p3p4 - ame**2*ampi**4*kp4*1D0 - 
     &         ampi**4*kp4*p1p2*1D0 - kp2*kp3*kp4*p1p3*1D0 - 
     &         ampi**2*kp2*kp4*p1p4*1D0 - 
     &         kp1*kp3*kp4*p2p3*1D0 - 
     &         ampi**2*kp4*p1p4*p2p3*1D0 - 
     &         ampi**2*kp1*kp4*p2p4*1D0 - 
     &         ampi**2*kp4*p1p3*p2p4*1D0 - 
     &         ame**2*kp3*kp4*p3p4*1D0 - 
     &         kp3*kp4*p1p2*p3p4*1D0 - 
     &         kp3*p1p3*p2p3*p3p4*1D0 - 
     &         kp3*p1p4*p2p4*p3p4*1D0 - 
     &         ame**2*kp3*p3p4**2*1D0 - 
     &         kp3*p1p2*p3p4**2*1D0 - kp3*kp4*p1p3*p2p3*2D0))/
     &     (kp3*kp4*(ame**2 + p1p2)**2)) + 
     &  (C1apx350pm1p*(ampi**2*kp1*kp2*kp4 + 
     &       ame**2*kp3**2*kp4 + ame**2*ampi**2*kp4**2 + 
     &       kp3**2*kp4*p1p2 + ampi**2*kp4**2*p1p2 + 
     &       kp3*kp4*p1p4*p2p3 + kp3*kp4*p1p3*p2p4 + 
     &       ampi**2*kp4*p1p4*p2p4 + kp1*kp2*kp3*p3p4 + 
     &       ame**2*kp3**2*p3p4 + kp3**2*p1p2*p3p4 + 
     &       kp3*p1p3*p2p3*p3p4 - 
     &       ame**2*ampi**2*kp3*kp4*1D0 - 
     &       ampi**2*kp3*kp4*p1p2*1D0 - 
     &       kp2*kp3*kp4*p1p3*1D0 - 
     &       ampi**2*kp2*kp4*p1p4*1D0 - 
     &       kp1*kp3*kp4*p2p3*1D0 - 
     &       ampi**2*kp4*p1p3*p2p3*1D0 - 
     &       ampi**2*kp1*kp4*p2p4*1D0 - 
     &       ame**2*kp3*kp4*p3p4*1D0 - 
     &       kp3*kp4*p1p2*p3p4*1D0 - kp2*kp3*p1p3*p3p4*1D0 - 
     &       kp1*kp3*p2p3*p3p4*1D0 - kp3*p1p4*p2p4*p3p4*1D0 + 
     &       kp3*kp4*p1p3*p2p3*2D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C0apx45x12pm1p*1D0*
     &     (ampi**2*kp1*kp2*kp3 + ame**2*ampi**2*kp3**2 + 
     &       ame**2*kp3*kp4**2 + ampi**2*kp3**2*p1p2 + 
     &       kp3*kp4**2*p1p2 + ampi**2*kp2*kp3*p1p4 + 
     &       ampi**2*kp3*p1p3*p2p3 + kp3*kp4*p1p4*p2p3 + 
     &       ampi**2*kp1*kp3*p2p4 + kp3*kp4*p1p3*p2p4 + 
     &       ampi**2*kp3*p1p4*p2p4 + 
     &       ame**2*ampi**2*kp3*p3p4 + 
     &       ame**2*ampi**2*kp4*p3p4 + kp1*kp2*kp4*p3p4 + 
     &       ampi**2*kp3*p1p2*p3p4 + ampi**2*kp4*p1p2*p3p4 + 
     &       kp4*p1p4*p2p3*p3p4 + kp4*p1p3*p2p4*p3p4 - 
     &       ame**2*ampi**4*kp3*1D0 - ampi**4*kp3*p1p2*1D0 - 
     &       ampi**2*kp2*kp3*p1p3*1D0 - 
     &       kp2*kp3*kp4*p1p4*1D0 - 
     &       ampi**2*kp1*kp3*p2p3*1D0 - 
     &       ampi**2*kp3*p1p4*p2p3*1D0 - 
     &       kp1*kp3*kp4*p2p4*1D0 - 
     &       ampi**2*kp3*p1p3*p2p4*1D0 - 
     &       ame**2*kp3*kp4*p3p4*1D0 - 
     &       kp3*kp4*p1p2*p3p4*1D0 - kp4*p1p3*p2p3*p3p4*1D0 - 
     &       kp4*p1p4*p2p4*p3p4*1D0 - 
     &       ame**2*kp4*p3p4**2*1D0 - kp4*p1p2*p3p4**2*1D0 - 
     &       kp3*kp4*p1p4*p2p4*2D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C1ax12x340ppp*(ampi**2*kp1*kp2*kp3 + 
     &       ampi**2*kp1*kp2*kp4 + ame**2*kp3**2*kp4 + 
     &       ame**2*kp3*kp4**2 + kp3**2*kp4*p1p2 + 
     &       kp3*kp4**2*p1p2 + ampi**2*kp2*kp4*p1p3 + 
     &       ampi**2*kp2*kp3*p1p4 + ampi**2*kp1*kp4*p2p3 + 
     &       ampi**2*kp4*p1p3*p2p3 + ampi**2*kp1*kp3*p2p4 + 
     &       ampi**2*kp3*p1p4*p2p4 + kp1*kp2*kp3*p3p4 + 
     &       kp1*kp2*kp4*p3p4 + kp2*kp4*p1p3*p3p4 + 
     &       kp2*kp3*p1p4*p3p4 + kp1*kp4*p2p3*p3p4 + 
     &       kp4*p1p3*p2p3*p3p4 + kp1*kp3*p2p4*p3p4 + 
     &       kp3*p1p4*p2p4*p3p4 - kp2*kp3*kp4*p1p3*1D0 - 
     &       kp2*kp3*kp4*p1p4*1D0 - kp1*kp3*kp4*p2p3*1D0 - 
     &       ampi**2*kp3*p1p3*p2p3*1D0 - 
     &       kp1*kp3*kp4*p2p4*1D0 - 
     &       ampi**2*kp4*p1p4*p2p4*1D0 - 
     &       kp3*p1p3*p2p3*p3p4*1D0 - 
     &       kp4*p1p4*p2p4*p3p4*1D0 - kp3*kp4*p1p3*p2p3*2D0 - 
     &       kp3*kp4*p1p4*p2p3*2D0 - kp3*kp4*p1p3*p2p4*2D0 - 
     &       kp3*kp4*p1p4*p2p4*2D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C11ax12x340ppp*2D0*
     &     (ampi**2*kp1*kp2*kp3 + ampi**2*kp1*kp2*kp4 + 
     &       ame**2*kp3**2*kp4 + ame**2*kp3*kp4**2 + 
     &       kp3**2*kp4*p1p2 + kp3*kp4**2*p1p2 + 
     &       ampi**2*kp2*kp4*p1p3 + ampi**2*kp2*kp3*p1p4 + 
     &       ampi**2*kp1*kp4*p2p3 + ampi**2*kp4*p1p3*p2p3 + 
     &       ampi**2*kp1*kp3*p2p4 + ampi**2*kp3*p1p4*p2p4 + 
     &       kp1*kp2*kp3*p3p4 + kp1*kp2*kp4*p3p4 + 
     &       kp2*kp4*p1p3*p3p4 + kp2*kp3*p1p4*p3p4 + 
     &       kp1*kp4*p2p3*p3p4 + kp4*p1p3*p2p3*p3p4 + 
     &       kp1*kp3*p2p4*p3p4 + kp3*p1p4*p2p4*p3p4 - 
     &       kp2*kp3*kp4*p1p3*1D0 - kp2*kp3*kp4*p1p4*1D0 - 
     &       kp1*kp3*kp4*p2p3*1D0 - 
     &       ampi**2*kp3*p1p3*p2p3*1D0 - 
     &       kp1*kp3*kp4*p2p4*1D0 - 
     &       ampi**2*kp4*p1p4*p2p4*1D0 - 
     &       kp3*p1p3*p2p3*p3p4*1D0 - 
     &       kp4*p1p4*p2p4*p3p4*1D0 - kp3*kp4*p1p3*p2p3*2D0 - 
     &       kp3*kp4*p1p4*p2p3*2D0 - kp3*kp4*p1p3*p2p4*2D0 - 
     &       kp3*kp4*p1p4*p2p4*2D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C1apx450pm1p*(ampi**2*kp1*kp2*kp3 + 
     &       ame**2*ampi**2*kp3**2 + ame**2*kp3*kp4**2 + 
     &       ampi**2*kp3**2*p1p2 + kp3*kp4**2*p1p2 + 
     &       ampi**2*kp3*p1p3*p2p3 + kp3*kp4*p1p4*p2p3 + 
     &       kp3*kp4*p1p3*p2p4 + kp1*kp2*kp4*p3p4 + 
     &       ame**2*kp4**2*p3p4 + kp4**2*p1p2*p3p4 + 
     &       kp4*p1p4*p2p4*p3p4 - 
     &       ame**2*ampi**2*kp3*kp4*1D0 - 
     &       ampi**2*kp3*kp4*p1p2*1D0 - 
     &       ampi**2*kp2*kp3*p1p3*1D0 - 
     &       kp2*kp3*kp4*p1p4*1D0 - 
     &       ampi**2*kp1*kp3*p2p3*1D0 - 
     &       kp1*kp3*kp4*p2p4*1D0 - 
     &       ampi**2*kp3*p1p4*p2p4*1D0 - 
     &       ame**2*kp3*kp4*p3p4*1D0 - 
     &       kp3*kp4*p1p2*p3p4*1D0 - kp2*kp4*p1p4*p3p4*1D0 - 
     &       kp4*p1p3*p2p3*p3p4*1D0 - kp1*kp4*p2p4*p3p4*1D0 + 
     &       kp3*kp4*p1p4*p2p4*2D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C12apx350pm1p*(ame**2*ampi**2*kp4**2 + 
     &       ampi**2*kp4**2*p1p2 + ampi**2*kp2*kp4*p1p3 + 
     &       ampi**2*kp1*kp4*p2p3 + ame**2*kp3**2*p3p4 + 
     &       kp3**2*p1p2*p3p4 + kp2*kp3*p1p4*p3p4 + 
     &       kp1*kp3*p2p4*p3p4 - ame**2*ampi**2*kp3*kp4*1D0 - 
     &       ampi**2*kp3*kp4*p1p2*1D0 - 
     &       ampi**2*kp2*kp4*p1p4*1D0 - 
     &       ampi**2*kp1*kp4*p2p4*1D0 - 
     &       ame**2*kp3*kp4*p3p4*1D0 - 
     &       kp3*kp4*p1p2*p3p4*1D0 - kp2*kp3*p1p3*p3p4*1D0 - 
     &       kp1*kp3*p2p3*p3p4*1D0 + 
     &       ampi**2*kp1*kp2*kp4*2D0 + 
     &       ame**2*kp3**2*kp4*2D0 + kp3**2*kp4*p1p2*2D0 - 
     &       kp2*kp3*kp4*p1p3*2D0 - kp1*kp3*kp4*p2p3*2D0 + 
     &       kp1*kp2*kp3*p3p4*2D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C12apx450pm1p*(ame**2*ampi**2*kp3**2 + 
     &       ampi**2*kp3**2*p1p2 + ampi**2*kp2*kp3*p1p4 + 
     &       ampi**2*kp1*kp3*p2p4 + ame**2*kp4**2*p3p4 + 
     &       kp4**2*p1p2*p3p4 + kp2*kp4*p1p3*p3p4 + 
     &       kp1*kp4*p2p3*p3p4 - ame**2*ampi**2*kp3*kp4*1D0 - 
     &       ampi**2*kp3*kp4*p1p2*1D0 - 
     &       ampi**2*kp2*kp3*p1p3*1D0 - 
     &       ampi**2*kp1*kp3*p2p3*1D0 - 
     &       ame**2*kp3*kp4*p3p4*1D0 - 
     &       kp3*kp4*p1p2*p3p4*1D0 - kp2*kp4*p1p4*p3p4*1D0 - 
     &       kp1*kp4*p2p4*p3p4*1D0 + 
     &       ampi**2*kp1*kp2*kp3*2D0 + 
     &       ame**2*kp3*kp4**2*2D0 + kp3*kp4**2*p1p2*2D0 - 
     &       kp2*kp3*kp4*p1p4*2D0 - kp1*kp3*kp4*p2p4*2D0 + 
     &       kp1*kp2*kp4*p3p4*2D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C1ax12x340ppp*1D0*
     &     (ampi**2*kp1*kp2*kp3 + ame**2*ampi**2*kp3**2 + 
     &       ampi**2*kp1*kp2*kp4 + ame**2*kp3**2*kp4 + 
     &       ame**2*ampi**2*kp4**2 + ame**2*kp3*kp4**2 + 
     &       ampi**2*kp3**2*p1p2 + kp3**2*kp4*p1p2 + 
     &       ampi**2*kp4**2*p1p2 + kp3*kp4**2*p1p2 + 
     &       ampi**2*kp3*p1p3*p2p3 + ampi**2*kp4*p1p4*p2p4 + 
     &       kp1*kp2*kp3*p3p4 + ame**2*kp3**2*p3p4 + 
     &       kp1*kp2*kp4*p3p4 + ame**2*kp4**2*p3p4 + 
     &       kp3**2*p1p2*p3p4 + kp4**2*p1p2*p3p4 + 
     &       kp3*p1p3*p2p3*p3p4 + kp4*p1p4*p2p4*p3p4 - 
     &       ampi**2*kp2*kp3*p1p3*1D0 - 
     &       kp2*kp3*kp4*p1p3*1D0 - 
     &       ampi**2*kp2*kp4*p1p4*1D0 - 
     &       kp2*kp3*kp4*p1p4*1D0 - 
     &       ampi**2*kp1*kp3*p2p3*1D0 - 
     &       kp1*kp3*kp4*p2p3*1D0 - 
     &       ampi**2*kp4*p1p3*p2p3*1D0 - 
     &       ampi**2*kp1*kp4*p2p4*1D0 - 
     &       kp1*kp3*kp4*p2p4*1D0 - 
     &       ampi**2*kp3*p1p4*p2p4*1D0 - 
     &       kp2*kp3*p1p3*p3p4*1D0 - kp2*kp4*p1p4*p3p4*1D0 - 
     &       kp1*kp3*p2p3*p3p4*1D0 - kp4*p1p3*p2p3*p3p4*1D0 - 
     &       kp1*kp4*p2p4*p3p4*1D0 - kp3*p1p4*p2p4*p3p4*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*2D0 - 
     &       ampi**2*kp3*kp4*p1p2*2D0 + 
     &       kp3*kp4*p1p3*p2p3*2D0 + kp3*kp4*p1p4*p2p3*2D0 + 
     &       kp3*kp4*p1p3*p2p4*2D0 + kp3*kp4*p1p4*p2p4*2D0 - 
     &       ame**2*kp3*kp4*p3p4*2D0 - kp3*kp4*p1p2*p3p4*2D0))
     &    /(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C12a0x34x12ppp*(ame**2*ampi**2*kp3**2 + 
     &       ame**2*ampi**2*kp4**2 + ampi**2*kp3**2*p1p2 + 
     &       ampi**2*kp4**2*p1p2 + ampi**2*kp2*kp4*p1p3 + 
     &       ampi**2*kp2*kp3*p1p4 + ampi**2*kp1*kp4*p2p3 + 
     &       ampi**2*kp1*kp3*p2p4 + ame**2*kp3**2*p3p4 + 
     &       ame**2*kp4**2*p3p4 + kp3**2*p1p2*p3p4 + 
     &       kp4**2*p1p2*p3p4 + kp2*kp4*p1p3*p3p4 + 
     &       kp2*kp3*p1p4*p3p4 + kp1*kp4*p2p3*p3p4 + 
     &       kp1*kp3*p2p4*p3p4 - ampi**2*kp2*kp3*p1p3*1D0 - 
     &       ampi**2*kp2*kp4*p1p4*1D0 - 
     &       ampi**2*kp1*kp3*p2p3*1D0 - 
     &       ampi**2*kp1*kp4*p2p4*1D0 - 
     &       kp2*kp3*p1p3*p3p4*1D0 - kp2*kp4*p1p4*p3p4*1D0 - 
     &       kp1*kp3*p2p3*p3p4*1D0 - kp1*kp4*p2p4*p3p4*1D0 + 
     &       ampi**2*kp1*kp2*kp3*2D0 + 
     &       ampi**2*kp1*kp2*kp4*2D0 - 
     &       ame**2*ampi**2*kp3*kp4*2D0 + 
     &       ame**2*kp3**2*kp4*2D0 + ame**2*kp3*kp4**2*2D0 - 
     &       ampi**2*kp3*kp4*p1p2*2D0 + kp3**2*kp4*p1p2*2D0 + 
     &       kp3*kp4**2*p1p2*2D0 - kp2*kp3*kp4*p1p3*2D0 - 
     &       kp2*kp3*kp4*p1p4*2D0 - kp1*kp3*kp4*p2p3*2D0 - 
     &       kp1*kp3*kp4*p2p4*2D0 + kp1*kp2*kp3*p3p4*2D0 + 
     &       kp1*kp2*kp4*p3p4*2D0 - ame**2*kp3*kp4*p3p4*2D0 - 
     &       kp3*kp4*p1p2*p3p4*2D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C12a0x12x34ppp*1D0*
     &     (ame**2*ampi**2*kp3**2 + ame**2*ampi**2*kp4**2 + 
     &       ampi**2*kp3**2*p1p2 + ampi**2*kp4**2*p1p2 + 
     &       ampi**2*kp2*kp4*p1p3 + ampi**2*kp2*kp3*p1p4 + 
     &       ampi**2*kp1*kp4*p2p3 + ampi**2*kp1*kp3*p2p4 + 
     &       ame**2*kp3**2*p3p4 + ame**2*kp4**2*p3p4 + 
     &       kp3**2*p1p2*p3p4 + kp4**2*p1p2*p3p4 + 
     &       kp2*kp4*p1p3*p3p4 + kp2*kp3*p1p4*p3p4 + 
     &       kp1*kp4*p2p3*p3p4 + kp1*kp3*p2p4*p3p4 - 
     &       ampi**2*kp2*kp3*p1p3*1D0 - 
     &       ampi**2*kp2*kp4*p1p4*1D0 - 
     &       ampi**2*kp1*kp3*p2p3*1D0 - 
     &       ampi**2*kp1*kp4*p2p4*1D0 - 
     &       kp2*kp3*p1p3*p3p4*1D0 - kp2*kp4*p1p4*p3p4*1D0 - 
     &       kp1*kp3*p2p3*p3p4*1D0 - kp1*kp4*p2p4*p3p4*1D0 + 
     &       ampi**2*kp1*kp2*kp3*2D0 + 
     &       ampi**2*kp1*kp2*kp4*2D0 - 
     &       ame**2*ampi**2*kp3*kp4*2D0 + 
     &       ame**2*kp3**2*kp4*2D0 + ame**2*kp3*kp4**2*2D0 - 
     &       ampi**2*kp3*kp4*p1p2*2D0 + kp3**2*kp4*p1p2*2D0 + 
     &       kp3*kp4**2*p1p2*2D0 - kp2*kp3*kp4*p1p3*2D0 - 
     &       kp2*kp3*kp4*p1p4*2D0 - kp1*kp3*kp4*p2p3*2D0 - 
     &       kp1*kp3*kp4*p2p4*2D0 + kp1*kp2*kp3*p3p4*2D0 + 
     &       kp1*kp2*kp4*p3p4*2D0 - ame**2*kp3*kp4*p3p4*2D0 - 
     &       kp3*kp4*p1p2*p3p4*2D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C12apx350pm1p*1D0*
     &     (kp3*p3p4*((ame**2 + p1p2)*(kp3 - kp4*1D0) + 
     &          kp2*(p1p4 - p1p3*1D0) + 
     &          kp1*(p2p4 - p2p3*1D0 + kp2*2D0)) + 
     &       kp4*(kp3**2*x12 + 
     &          ampi**2*(kp4*(ame**2 + p1p2) + 
     &             kp2*(p1p3 - p1p4*1D0)) - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) + 
     &             kp2*p1p3*2D0) + 
     &          kp1*(-(ampi**2*p2p4*1D0) + ampi**2*kp2*2D0 + 
     &             p2p3*(ampi**2 - kp3*2D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (A0ap*2.5D-1*(-(kp4*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(kp3*(ame**2 + p1p2)*(ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(ame**2*kp4**2 + kp4**2*p1p2 + 
     &          ampi**2*kp3*(ame**2 + p1p2) + kp4*p1p4*p2p3 - 
     &          kp2*kp4*p1p4*1D0 - ampi**2*p1p4*p2p3*1D0 - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0) + 
     &             ampi**2*(kp2 - p2p3*1D0)) + 
     &          p2p4*(p1p3*(kp4 - ampi**2*1D0) + 
     &             p1p4*(ampi**2 - kp4*2D0)))))/
     &   (kp3*kp4**3*(ame**2 + p1p2)**2) + 
     &  (A0am1*2.5D-1*(-(kp4*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(kp3*(ame**2 + p1p2)*(ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(ame**2*kp4**2 + kp4**2*p1p2 + 
     &          ampi**2*kp3*(ame**2 + p1p2) + kp4*p1p4*p2p3 - 
     &          kp2*kp4*p1p4*1D0 - ampi**2*p1p4*p2p3*1D0 - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0) + 
     &             ampi**2*(kp2 - p2p3*1D0)) + 
     &          p2p4*(p1p3*(kp4 - ampi**2*1D0) + 
     &             p1p4*(ampi**2 - kp4*2D0)))))/
     &   (ampi**2*kp3*kp4**2*(ame**2 + p1p2)**2) + 
     &  (A0am1*2.5D-1*(-(kp4*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(kp3*(ame**2 + p1p2)*(ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(ame**2*kp4**2 + kp4**2*p1p2 + 
     &          ampi**2*kp3*(ame**2 + p1p2) + kp4*p1p4*p2p3 - 
     &          kp2*kp4*p1p4*1D0 - ampi**2*p1p4*p2p3*1D0 - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0) + 
     &             ampi**2*(kp2 - p2p3*1D0)) + 
     &          p2p4*(p1p3*(kp4 - ampi**2*1D0) + 
     &             p1p4*(ampi**2 - kp4*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2*x45) - 
     &  (A0ap*2.5D-1*(-(kp4*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(kp3*(ame**2 + p1p2)*(ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(ame**2*kp4**2 + kp4**2*p1p2 + 
     &          ampi**2*kp3*(ame**2 + p1p2) + kp4*p1p4*p2p3 - 
     &          kp2*kp4*p1p4*1D0 - ampi**2*p1p4*p2p3*1D0 - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0) + 
     &             ampi**2*(kp2 - p2p3*1D0)) + 
     &          p2p4*(p1p3*(kp4 - ampi**2*1D0) + 
     &             p1p4*(ampi**2 - kp4*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2*x45) - 
     &  (A0am1*5.D-1*(-(kp4*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(kp3*(ame**2 + p1p2)*(ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(ame**2*kp4**2 + kp4**2*p1p2 + 
     &          ampi**2*kp3*(ame**2 + p1p2) + kp4*p1p4*p2p3 - 
     &          kp2*kp4*p1p4*1D0 - ampi**2*p1p4*p2p3*1D0 - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0) + 
     &             ampi**2*(kp2 - p2p3*1D0)) + 
     &          p2p4*(p1p3*(kp4 - ampi**2*1D0) + 
     &             p1p4*(ampi**2 - kp4*2D0)))))/
     &   (kp3*kp4**3*(ame**2 + p1p2)**2) - 
     &  (B0ax45pm1*2.5D-1*
     &     (-(m1**2*1D0) + ampi**2*2D0 + kp4*2D0)*
     &     (-(kp4*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(kp3*(ame**2 + p1p2)*(ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(ame**2*kp4**2 + kp4**2*p1p2 + 
     &          ampi**2*kp3*(ame**2 + p1p2) + kp4*p1p4*p2p3 - 
     &          kp2*kp4*p1p4*1D0 - ampi**2*p1p4*p2p3*1D0 - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0) + 
     &             ampi**2*(kp2 - p2p3*1D0)) + 
     &          p2p4*(p1p3*(kp4 - ampi**2*1D0) + 
     &             p1p4*(ampi**2 - kp4*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2*x45) - 
     &  (C12apx450pm1p*1D0*
     &     (kp4*p3p4*(-((ame**2 + p1p2)*1D0*
     &             (kp3 - kp4*1D0)) + kp2*(p1p3 - p1p4*1D0) + 
     &          kp1*(p2p3 - p2p4*1D0 + kp2*2D0)) + 
     &       kp3*(ampi**2*kp3*(ame**2 + p1p2) - 
     &          ame**2*ampi**2*kp4*1D0 - 
     &          ampi**2*kp4*p1p2*1D0 - ampi**2*kp2*p1p3*1D0 + 
     &          ame**2*kp4**2*2D0 + kp4**2*p1p2*2D0 + 
     &          kp2*p1p4*(ampi**2 - kp4*2D0) + 
     &          kp1*(ampi**2*(-(p2p3*1D0) + kp2*2D0) + 
     &             p2p4*(ampi**2 - kp4*2D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (A0ap*2.5D-1*(kp4*(ame**2 + p1p2)*p3p4**2 - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(kp2*kp4*p1p4 + ampi**2*p1p4*p2p3 - 
     &          ame**2*kp4**2*1D0 - kp4**2*p1p2*1D0 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          kp4*p1p4*p2p3*1D0 + 
     &          kp1*(p2p4*(kp4 - ampi**2*1D0) + 
     &             ampi**2*(p2p3 - kp2*1D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          p2p4*(p1p3*(ampi**2 - kp4*1D0) - 
     &             p1p4*1D0*(ampi**2 - kp4*2D0)))))/
     &   (ampi**2*kp3*kp4**2*(ame**2 + p1p2)**2) + 
     &  (A0am1*(ampi**2 + kp4)*5.D-1*
     &     (kp4*(ame**2 + p1p2)*p3p4**2 - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(kp2*kp4*p1p4 + ampi**2*p1p4*p2p3 - 
     &          ame**2*kp4**2*1D0 - kp4**2*p1p2*1D0 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          kp4*p1p4*p2p3*1D0 + 
     &          kp1*(p2p4*(kp4 - ampi**2*1D0) + 
     &             ampi**2*(p2p3 - kp2*1D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          p2p4*(p1p3*(ampi**2 - kp4*1D0) - 
     &             p1p4*1D0*(ampi**2 - kp4*2D0)))))/
     &   (ampi**2*kp3*kp4**2*(ame**2 + p1p2)**2*x45) - 
     &  (A0ap*(ampi**2 + kp4)*5.D-1*
     &     (kp4*(ame**2 + p1p2)*p3p4**2 - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(kp2*kp4*p1p4 + ampi**2*p1p4*p2p3 - 
     &          ame**2*kp4**2*1D0 - kp4**2*p1p2*1D0 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          kp4*p1p4*p2p3*1D0 + 
     &          kp1*(p2p4*(kp4 - ampi**2*1D0) + 
     &             ampi**2*(p2p3 - kp2*1D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          p2p4*(p1p3*(ampi**2 - kp4*1D0) - 
     &             p1p4*1D0*(ampi**2 - kp4*2D0)))))/
     &   (ampi**2*kp3*kp4**2*(ame**2 + p1p2)**2*x45) + 
     &  (B0appm1*2.5D-1*(-(m1**2*1D0) + ampi**2*2D0)*
     &     (kp4*(ame**2 + p1p2)*p3p4**2 - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(kp2*kp4*p1p4 + ampi**2*p1p4*p2p3 - 
     &          ame**2*kp4**2*1D0 - kp4**2*p1p2*1D0 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          kp4*p1p4*p2p3*1D0 + 
     &          kp1*(p2p4*(kp4 - ampi**2*1D0) + 
     &             ampi**2*(p2p3 - kp2*1D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          p2p4*(p1p3*(ampi**2 - kp4*1D0) - 
     &             p1p4*1D0*(ampi**2 - kp4*2D0)))))/
     &   (ampi**2*kp3*kp4**2*(ame**2 + p1p2)**2) - 
     &  (A0am1*2.5D-1*(kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)))))/
     &   (ampi**2*kp3*kp4**2*(ame**2 + p1p2)**2) + 
     &  (A0ap*2.5D-1*(kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)))))/
     &   (ampi**2*kp3*kp4**2*(ame**2 + p1p2)**2) + 
     &  (A0am1*2.5D-1*(kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2*x45) - 
     &  (A0ap*2.5D-1*(kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2*x45) - 
     &  (C11apx12x45m1pp*(kp3 + p3p4)*1D0*
     &     (kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2) + 
     &  (C11ax45x12pm1pp*
     &     (kp1 + kp2 + p1p4 + p2p4 - kp3*1D0 - p3p4*1D0)*
     &     (kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2) - 
     &  (C0apx12x45m1pp*(kp3 + p3p4)*2D0*
     &     (kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2) + 
     &  (B0ax45pm1*2.5D-1*
     &     (-(m1**2*1D0) + (ame**2 + ampi**2)*2D0 + 
     &       p1p2*2D0 - p1p3*2D0 - p2p3*2D0)*
     &     (kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2*x45) + 
     &  (C12apx12x45m1pp*
     &     (kp1 + kp2 + p1p4 + p2p4 - kp3*2D0 - p3p4*2D0)*
     &     (kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2) + 
     &  (A0ap*2.5D-1*(kp3*(ame**2 + p1p2)*p3p4**2 - 
     &       p3p4*1D0*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(-(kp3**2*(ame**2 + p1p2)*1D0) + 
     &          kp1*(ampi**2*p2p4 - ampi**2*kp2*1D0 + 
     &             p2p3*(kp3 - ampi**2*1D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*
     &              (kp2 + p2p3 - p2p4*1D0)) + 
     &          kp3*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &             p2p3*(-(p1p4*1D0) + p1p3*2D0)))))/
     &   (ampi**2*kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (A0am1*(ampi**2 + kp3)*5.D-1*
     &     (kp3*(ame**2 + p1p2)*p3p4**2 - 
     &       p3p4*1D0*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(-(kp3**2*(ame**2 + p1p2)*1D0) + 
     &          kp1*(ampi**2*p2p4 - ampi**2*kp2*1D0 + 
     &             p2p3*(kp3 - ampi**2*1D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*
     &              (kp2 + p2p3 - p2p4*1D0)) + 
     &          kp3*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &             p2p3*(-(p1p4*1D0) + p1p3*2D0)))))/
     &   (ampi**2*kp3**2*kp4*(ame**2 + p1p2)**2*x35) - 
     &  (A0ap*(ampi**2 + kp3)*5.D-1*
     &     (kp3*(ame**2 + p1p2)*p3p4**2 - 
     &       p3p4*1D0*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(-(kp3**2*(ame**2 + p1p2)*1D0) + 
     &          kp1*(ampi**2*p2p4 - ampi**2*kp2*1D0 + 
     &             p2p3*(kp3 - ampi**2*1D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*
     &              (kp2 + p2p3 - p2p4*1D0)) + 
     &          kp3*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &             p2p3*(-(p1p4*1D0) + p1p3*2D0)))))/
     &   (ampi**2*kp3**2*kp4*(ame**2 + p1p2)**2*x35) + 
     &  (B0appm1*2.5D-1*(-(m1**2*1D0) + ampi**2*2D0)*
     &     (kp3*(ame**2 + p1p2)*p3p4**2 - 
     &       p3p4*1D0*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(-(kp3**2*(ame**2 + p1p2)*1D0) + 
     &          kp1*(ampi**2*p2p4 - ampi**2*kp2*1D0 + 
     &             p2p3*(kp3 - ampi**2*1D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*
     &              (kp2 + p2p3 - p2p4*1D0)) + 
     &          kp3*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &             p2p3*(-(p1p4*1D0) + p1p3*2D0)))))/
     &   (ampi**2*kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (A0am1*2.5D-1*(-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(p1p4*(kp2 + p2p3) + p1p3*p2p4 + 
     &             (ame**2 + p1p2)*(ampi**2 - kp4*1D0) - 
     &             p1p3*p2p3*2D0 + kp1*(p2p4 + kp2*2D0))) + 
     &       kp4*(kp3**2*x12 - 
     &          kp3*p1p3*2D0*(kp2 + p2p3*2D0) + 
     &          kp1*(-(ampi**2*p2p4*1D0) + ampi**2*kp2*2D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)*2D0) - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + p1p3*p2p4 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             (kp2 + p2p3)*1D0*(-(p1p4*1D0) + p1p3*2D0)))
     &       ))/(kp3**2*kp4*(ame**2 + p1p2)**2*x35) - 
     &  (A0ap*2.5D-1*(-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(p1p4*(kp2 + p2p3) + p1p3*p2p4 + 
     &             (ame**2 + p1p2)*(ampi**2 - kp4*1D0) - 
     &             p1p3*p2p3*2D0 + kp1*(p2p4 + kp2*2D0))) + 
     &       kp4*(kp3**2*x12 - 
     &          kp3*p1p3*2D0*(kp2 + p2p3*2D0) + 
     &          kp1*(-(ampi**2*p2p4*1D0) + ampi**2*kp2*2D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)*2D0) - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + p1p3*p2p4 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             (kp2 + p2p3)*1D0*(-(p1p4*1D0) + p1p3*2D0)))
     &       ))/(kp3**2*kp4*(ame**2 + p1p2)**2*x35) - 
     &  (C00apx12x35m1pp*1D0*
     &     (-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(p1p4*(kp2 + p2p3) + p1p3*p2p4 + 
     &             (ame**2 + p1p2)*(ampi**2 - kp4*1D0) - 
     &             p1p3*p2p3*2D0 + kp1*(p2p4 + kp2*2D0))) + 
     &       kp4*(kp3**2*x12 - 
     &          kp3*p1p3*2D0*(kp2 + p2p3*2D0) + 
     &          kp1*(-(ampi**2*p2p4*1D0) + ampi**2*kp2*2D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)*2D0) - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + p1p3*p2p4 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             (kp2 + p2p3)*1D0*(-(p1p4*1D0) + p1p3*2D0)))
     &       ))/(kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (A0ap*2.5D-1*(-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(kp3**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*kp2 - ampi**2*p2p4*1D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)) + 
     &          ampi**2*(kp4*(ame**2 + p1p2) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &             (p1p3 - p1p4*1D0)*(kp2 + p2p3 - p2p4*1D0))
     &           + kp3*(p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             p1p3*(p2p4 - p2p3*2D0)))))/
     &   (kp3**3*kp4*(ame**2 + p1p2)**2) + 
     &  (A0am1*2.5D-1*(-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(kp3**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*kp2 - ampi**2*p2p4*1D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)) + 
     &          ampi**2*(kp4*(ame**2 + p1p2) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &             (p1p3 - p1p4*1D0)*(kp2 + p2p3 - p2p4*1D0))
     &           + kp3*(p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             p1p3*(p2p4 - p2p3*2D0)))))/
     &   (ampi**2*kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (A0am1*2.5D-1*(-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(kp3**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*kp2 - ampi**2*p2p4*1D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)) + 
     &          ampi**2*(kp4*(ame**2 + p1p2) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &             (p1p3 - p1p4*1D0)*(kp2 + p2p3 - p2p4*1D0))
     &           + kp3*(p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             p1p3*(p2p4 - p2p3*2D0)))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2*x35) - 
     &  (A0ap*2.5D-1*(-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(kp3**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*kp2 - ampi**2*p2p4*1D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)) + 
     &          ampi**2*(kp4*(ame**2 + p1p2) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &             (p1p3 - p1p4*1D0)*(kp2 + p2p3 - p2p4*1D0))
     &           + kp3*(p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             p1p3*(p2p4 - p2p3*2D0)))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2*x35) - 
     &  (A0am1*5.D-1*(-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(kp3**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*kp2 - ampi**2*p2p4*1D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)) + 
     &          ampi**2*(kp4*(ame**2 + p1p2) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &             (p1p3 - p1p4*1D0)*(kp2 + p2p3 - p2p4*1D0))
     &           + kp3*(p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             p1p3*(p2p4 - p2p3*2D0)))))/
     &   (kp3**3*kp4*(ame**2 + p1p2)**2) - 
     &  (B0ax35pm1*2.5D-1*
     &     (-(m1**2*1D0) + ampi**2*2D0 + kp3*2D0)*
     &     (-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(kp3**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*kp2 - ampi**2*p2p4*1D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)) + 
     &          ampi**2*(kp4*(ame**2 + p1p2) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &             (p1p3 - p1p4*1D0)*(kp2 + p2p3 - p2p4*1D0))
     &           + kp3*(p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             p1p3*(p2p4 - p2p3*2D0)))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2*x35) - 
     &  (A0am1*2.5D-1*(kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)))))/
     &   (ampi**2*kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (A0ap*2.5D-1*(kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)))))/
     &   (ampi**2*kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (A0am1*2.5D-1*(kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2*x35) - 
     &  (A0ap*2.5D-1*(kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2*x35) - 
     &  (C11apx12x35m1pp*(kp4 + p3p4)*1D0*
     &     (kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (C11ax35x12pm1pp*
     &     (kp1 + kp2 + p1p3 + p2p3 - kp4*1D0 - p3p4*1D0)*
     &     (kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2) - 
     &  (C0apx12x35m1pp*(kp4 + p3p4)*2D0*
     &     (kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (B0ax35pm1*2.5D-1*
     &     (-(m1**2*1D0) + (ame**2 + ampi**2)*2D0 + 
     &       p1p2*2D0 - p1p4*2D0 - p2p4*2D0)*
     &     (kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2*x35) + 
     &  (C12apx12x35m1pp*
     &     (kp1 + kp2 + p1p3 + p2p3 - kp4*2D0 - p3p4*2D0)*
     &     (kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (A0am1*2.5D-1*(ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*p3p4*(kp2*p1p3 + p1p4*p2p3 + 
     &          (ame**2 + p1p2)*(ampi**2 - p3p4*1D0) + 
     &          kp1*(p2p3 + kp2*2D0) + p2p4*(p1p3 - p1p4*2D0))
     &         + kp3*(ampi**2*(ame**2 + p1p2)*p3p4 + 
     &          kp4**2*x12 + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0)*2D0 + 
     &             ampi**2*(-(p2p3*1D0) + kp2*2D0)) - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + p1p3*p2p4 + 
     &             kp2*(p1p3 - p1p4*2D0) + 
     &             p1p4*(p2p3 - p2p4*2D0)) + 
     &          kp4*(-((ame**2 + p1p2)*p3p4*1D0) - 
     &             p1p4*2D0*(kp2 + p2p4*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2*x45) - 
     &  (A0ap*2.5D-1*(ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*p3p4*(kp2*p1p3 + p1p4*p2p3 + 
     &          (ame**2 + p1p2)*(ampi**2 - p3p4*1D0) + 
     &          kp1*(p2p3 + kp2*2D0) + p2p4*(p1p3 - p1p4*2D0))
     &         + kp3*(ampi**2*(ame**2 + p1p2)*p3p4 + 
     &          kp4**2*x12 + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0)*2D0 + 
     &             ampi**2*(-(p2p3*1D0) + kp2*2D0)) - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + p1p3*p2p4 + 
     &             kp2*(p1p3 - p1p4*2D0) + 
     &             p1p4*(p2p3 - p2p4*2D0)) + 
     &          kp4*(-((ame**2 + p1p2)*p3p4*1D0) - 
     &             p1p4*2D0*(kp2 + p2p4*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2*x45) - 
     &  (C00apx12x45m1pp*1D0*
     &     (ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*p3p4*(kp2*p1p3 + p1p4*p2p3 + 
     &          (ame**2 + p1p2)*(ampi**2 - p3p4*1D0) + 
     &          kp1*(p2p3 + kp2*2D0) + p2p4*(p1p3 - p1p4*2D0))
     &         + kp3*(ampi**2*(ame**2 + p1p2)*p3p4 + 
     &          kp4**2*x12 + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0)*2D0 + 
     &             ampi**2*(-(p2p3*1D0) + kp2*2D0)) - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + p1p3*p2p4 + 
     &             kp2*(p1p3 - p1p4*2D0) + 
     &             p1p4*(p2p3 - p2p4*2D0)) + 
     &          kp4*(-((ame**2 + p1p2)*p3p4*1D0) - 
     &             p1p4*2D0*(kp2 + p2p4*2D0)))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2) - 
     &  (C1apx12x45m1pp*(kp3 + p3p4)*
     &     (kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0))))*3D0)/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2) - 
     &  (C1apx12x35m1pp*(kp4 + p3p4)*
     &     (kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0))))*3D0)/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (C2ax45px12pm1p*(ame**2*ampi**4*kp3 + 
     &       ame**2*ampi**4*kp4 + ampi**2*kp1*kp2*kp4 + 
     &       ame**2*kp3**2*kp4 + ampi**4*kp3*p1p2 + 
     &       ampi**4*kp4*p1p2 + kp3**2*kp4*p1p2 + 
     &       ampi**2*kp2*kp3*p1p3 + ampi**2*kp1*kp3*p2p3 + 
     &       ampi**2*kp3*p1p4*p2p3 + ampi**2*kp4*p1p4*p2p3 + 
     &       ampi**2*kp3*p1p3*p2p4 + ampi**2*kp4*p1p3*p2p4 + 
     &       kp1*kp2*kp3*p3p4 + ame**2*kp3**2*p3p4 + 
     &       ame**2*kp3*kp4*p3p4 + kp3**2*p1p2*p3p4 + 
     &       kp3*kp4*p1p2*p3p4 + kp2*kp4*p1p3*p3p4 + 
     &       kp2*kp3*p1p4*p3p4 + kp1*kp4*p2p3*p3p4 + 
     &       kp3*p1p3*p2p3*p3p4 + kp1*kp3*p2p4*p3p4 + 
     &       kp3*p1p4*p2p4*p3p4 + ame**2*kp3*p3p4**2 + 
     &       ame**2*kp4*p3p4**2 + kp3*p1p2*p3p4**2 + 
     &       kp4*p1p2*p3p4**2 - ame**2*ampi**2*kp3**2*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*1D0 - 
     &       ampi**2*kp3**2*p1p2*1D0 - 
     &       ampi**2*kp3*kp4*p1p2*1D0 - 
     &       kp2*kp3*kp4*p1p3*1D0 - kp1*kp3*kp4*p2p3*1D0 - 
     &       ampi**2*kp4*p1p3*p2p3*1D0 - 
     &       ampi**2*kp4*p1p4*p2p4*1D0 - 
     &       kp2*kp3*p1p3*p3p4*1D0 - kp1*kp3*p2p3*p3p4*1D0 - 
     &       kp3*p1p4*p2p3*p3p4*1D0 - 
     &       kp4*p1p4*p2p3*p3p4*1D0 - 
     &       kp3*p1p3*p2p4*p3p4*1D0 - 
     &       kp4*p1p3*p2p4*p3p4*1D0 - 
     &       ampi**2*kp3*p1p3*p2p3*2D0 + 
     &       kp3*kp4*p1p3*p2p3*2D0 - 
     &       ame**2*ampi**2*kp3*p3p4*2D0 - 
     &       ame**2*ampi**2*kp4*p3p4*2D0 - 
     &       ampi**2*kp3*p1p2*p3p4*2D0 - 
     &       ampi**2*kp4*p1p2*p3p4*2D0 + 
     &       kp4*p1p3*p2p3*p3p4*2D0 - kp3*kp4*p1p4*p2p3*3D0 - 
     &       kp3*kp4*p1p3*p2p4*3D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C2ax35px12pm1p*(ame**2*ampi**4*kp3 + 
     &       ampi**2*kp1*kp2*kp3 + ame**2*ampi**4*kp4 + 
     &       ame**2*kp3*kp4**2 + ampi**4*kp3*p1p2 + 
     &       ampi**4*kp4*p1p2 + kp3*kp4**2*p1p2 + 
     &       ampi**2*kp2*kp4*p1p4 + ampi**2*kp3*p1p4*p2p3 + 
     &       ampi**2*kp4*p1p4*p2p3 + ampi**2*kp1*kp4*p2p4 + 
     &       ampi**2*kp3*p1p3*p2p4 + ampi**2*kp4*p1p3*p2p4 + 
     &       kp1*kp2*kp4*p3p4 + ame**2*kp3*kp4*p3p4 + 
     &       ame**2*kp4**2*p3p4 + kp3*kp4*p1p2*p3p4 + 
     &       kp4**2*p1p2*p3p4 + kp2*kp4*p1p3*p3p4 + 
     &       kp2*kp3*p1p4*p3p4 + kp1*kp4*p2p3*p3p4 + 
     &       kp4*p1p3*p2p3*p3p4 + kp1*kp3*p2p4*p3p4 + 
     &       kp4*p1p4*p2p4*p3p4 + ame**2*kp3*p3p4**2 + 
     &       ame**2*kp4*p3p4**2 + kp3*p1p2*p3p4**2 + 
     &       kp4*p1p2*p3p4**2 - ame**2*ampi**2*kp3*kp4*1D0 - 
     &       ame**2*ampi**2*kp4**2*1D0 - 
     &       ampi**2*kp3*kp4*p1p2*1D0 - 
     &       ampi**2*kp4**2*p1p2*1D0 - kp2*kp3*kp4*p1p4*1D0 - 
     &       ampi**2*kp3*p1p3*p2p3*1D0 - 
     &       kp1*kp3*kp4*p2p4*1D0 - 
     &       ampi**2*kp3*p1p4*p2p4*1D0 - 
     &       kp2*kp4*p1p4*p3p4*1D0 - kp3*p1p4*p2p3*p3p4*1D0 - 
     &       kp4*p1p4*p2p3*p3p4*1D0 - kp1*kp4*p2p4*p3p4*1D0 - 
     &       kp3*p1p3*p2p4*p3p4*1D0 - 
     &       kp4*p1p3*p2p4*p3p4*1D0 - 
     &       ampi**2*kp4*p1p4*p2p4*2D0 + 
     &       kp3*kp4*p1p4*p2p4*2D0 - 
     &       ame**2*ampi**2*kp3*p3p4*2D0 - 
     &       ame**2*ampi**2*kp4*p3p4*2D0 - 
     &       ampi**2*kp3*p1p2*p3p4*2D0 - 
     &       ampi**2*kp4*p1p2*p3p4*2D0 + 
     &       kp3*p1p4*p2p4*p3p4*2D0 - kp3*kp4*p1p4*p2p3*3D0 - 
     &       kp3*kp4*p1p3*p2p4*3D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C1ax45x12pm1pp*(kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0))))*
     &     (kp1 + kp2 + p1p4 + p2p4 - kp3*3D0 - p3p4*3D0))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2) + 
     &  (C1ax35x12pm1pp*(kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0))))*
     &     (kp1 + kp2 + p1p3 + p2p3 - kp4*3D0 - p3p4*3D0))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (B0appm1*(kp4*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C00apx12x35m1pp*
     &     (kp4*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C00apx12x45m1pp*
     &     (kp4*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C00apx350pm1p*(kp4*
     &        (kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C00apx450pm1p*(kp4*
     &        (kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (B0ax34pp*5.D-1*(kp4*
     &        (kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (B0ax35pm1*1D0*(kp4*
     &        (kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (B0ax45pm1*1D0*(kp4*
     &        (kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (B0appm1*2.5D-1*(kp4*(ame**2 + p1p2)*p3p4**2 - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(kp2*kp4*p1p4 + ampi**2*p1p4*p2p3 - 
     &          ame**2*kp4**2*1D0 - kp4**2*p1p2*1D0 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          kp4*p1p4*p2p3*1D0 + 
     &          kp1*(p2p4*(kp4 - ampi**2*1D0) + 
     &             ampi**2*(p2p3 - kp2*1D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          p2p4*(p1p3*(ampi**2 - kp4*1D0) - 
     &             p1p4*1D0*(ampi**2 - kp4*2D0))))*
     &     (-(m1**2*1D0) + ampi**2*4D0))/
     &   (ampi**2*kp3*kp4**2*(ame**2 + p1p2)**2) - 
     &  (B0appm1*2.5D-1*(kp4*(ame**2 + p1p2)*p3p4**2 + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 - 
     &          ampi**2*kp3*(ame**2 + p1p2)*1D0 - 
     &          ampi**2*p1p3*p2p3*2D0 - kp4*p1p4*p2p3*2D0 + 
     &          p1p3*p2p4*(ampi**2 - kp4*2D0)) - 
     &       p3p4*1D0*(kp3*(ame**2 + p1p2)*
     &           (ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0))))*
     &     (-(m1**2*1D0) + ampi**2*4D0))/
     &   (ampi**2*kp3*kp4**2*(ame**2 + p1p2)**2) - 
     &  (B0appm1*2.5D-1*(kp3*(ame**2 + p1p2)*p3p4**2 - 
     &       p3p4*1D0*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(-(kp3**2*(ame**2 + p1p2)*1D0) + 
     &          kp1*(ampi**2*p2p4 - ampi**2*kp2*1D0 + 
     &             p2p3*(kp3 - ampi**2*1D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*
     &              (kp2 + p2p3 - p2p4*1D0)) + 
     &          kp3*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &             p2p3*(-(p1p4*1D0) + p1p3*2D0))))*
     &     (-(m1**2*1D0) + ampi**2*4D0))/
     &   (ampi**2*kp3**2*kp4*(ame**2 + p1p2)**2) - 
     &  (B0appm1*2.5D-1*(kp3*(ame**2 + p1p2)*p3p4**2 + 
     &       kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ampi**2*kp4*(ame**2 + p1p2)*1D0 + 
     &          p2p4*(ampi**2*kp1 - ampi**2*p1p4*2D0 + 
     &             p1p3*(ampi**2 - kp3*2D0)) + 
     &          p1p4*(ampi**2*kp2 + p2p3*(ampi**2 - kp3*2D0)))
     &         - p3p4*1D0*
     &        (ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p2p4*(p1p3 - kp1*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0))))*
     &     (-(m1**2*1D0) + ampi**2*4D0))/
     &   (ampi**2*kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (C1apx350pm1p*5.D-1*
     &     (-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(kp3**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*kp2 - ampi**2*p2p4*1D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)) + 
     &          ampi**2*(kp4*(ame**2 + p1p2) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &             (p1p3 - p1p4*1D0)*(kp2 + p2p3 - p2p4*1D0))
     &           + kp3*(p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             p1p3*(p2p4 - p2p3*2D0))))*
     &     (-(m1**2*1D0) + kp3*2D0 + ampi**2*4D0))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2) + 
     &  (C1apx450pm1p*5.D-1*
     &     (-(kp4*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(kp3*(ame**2 + p1p2)*(ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(ame**2*kp4**2 + kp4**2*p1p2 + 
     &          ampi**2*kp3*(ame**2 + p1p2) + kp4*p1p4*p2p3 - 
     &          kp2*kp4*p1p4*1D0 - ampi**2*p1p4*p2p3*1D0 - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0) + 
     &             ampi**2*(kp2 - p2p3*1D0)) + 
     &          p2p4*(p1p3*(kp4 - ampi**2*1D0) + 
     &             p1p4*(ampi**2 - kp4*2D0))))*
     &     (-(m1**2*1D0) + kp4*2D0 + ampi**2*4D0))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2) - 
     &  (B0ax35pm1*2.5D-1*
     &     (-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(kp3**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*kp2 - ampi**2*p2p4*1D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)) + 
     &          ampi**2*(kp4*(ame**2 + p1p2) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &             (p1p3 - p1p4*1D0)*(kp2 + p2p3 - p2p4*1D0))
     &           + kp3*(p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             p1p3*(p2p4 - p2p3*2D0))))*
     &     (-(m1**2*1D0) + ampi**2*4D0 + kp3*4D0))/
     &   (kp3**3*kp4*(ame**2 + p1p2)**2) - 
     &  (B0ax45pm1*2.5D-1*
     &     (-(kp4*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(kp3*(ame**2 + p1p2)*(ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(ame**2*kp4**2 + kp4**2*p1p2 + 
     &          ampi**2*kp3*(ame**2 + p1p2) + kp4*p1p4*p2p3 - 
     &          kp2*kp4*p1p4*1D0 - ampi**2*p1p4*p2p3*1D0 - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0) + 
     &             ampi**2*(kp2 - p2p3*1D0)) + 
     &          p2p4*(p1p3*(kp4 - ampi**2*1D0) + 
     &             p1p4*(ampi**2 - kp4*2D0))))*
     &     (-(m1**2*1D0) + ampi**2*4D0 + kp4*4D0))/
     &   (kp3*kp4**3*(ame**2 + p1p2)**2) - 
     &  (C22apx35x12pm1p*1D0*
     &     (ame**2*ampi**2*kp4**2 + ampi**2*kp4**2*p1p2 + 
     &       ame**2*ampi**2*kp3*p3p4 + 
     &       ame**2*ampi**2*kp4*p3p4 + 
     &       ampi**2*kp3*p1p2*p3p4 + ampi**2*kp4*p1p2*p3p4 + 
     &       kp2*kp3*p1p4*p3p4 + kp3*p1p4*p2p3*p3p4 + 
     &       kp1*kp3*p2p4*p3p4 + kp3*p1p3*p2p4*p3p4 - 
     &       ame**2*ampi**4*kp4*1D0 - ampi**4*kp4*p1p2*1D0 - 
     &       ampi**2*kp2*kp4*p1p4*1D0 - 
     &       ampi**2*kp4*p1p4*p2p3*1D0 - 
     &       ampi**2*kp1*kp4*p2p4*1D0 - 
     &       ampi**2*kp4*p1p3*p2p4*1D0 - 
     &       ame**2*kp3*kp4*p3p4*1D0 - 
     &       kp3*kp4*p1p2*p3p4*1D0 - ame**2*kp3*p3p4**2*1D0 - 
     &       kp3*p1p2*p3p4**2*1D0 + ampi**2*kp1*kp2*kp4*2D0 + 
     &       ame**2*kp3**2*kp4*2D0 + kp3**2*kp4*p1p2*2D0 + 
     &       ampi**2*kp2*kp4*p1p3*2D0 - 
     &       kp2*kp3*kp4*p1p3*2D0 + 
     &       ampi**2*kp1*kp4*p2p3*2D0 - 
     &       kp1*kp3*kp4*p2p3*2D0 + 
     &       ampi**2*kp4*p1p3*p2p3*2D0 + 
     &       kp1*kp2*kp3*p3p4*2D0 - kp3*p1p3*p2p3*p3p4*2D0 - 
     &       kp3*kp4*p1p3*p2p3*4D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C11apx350pm1p*(ame**2*ampi**4*kp4 + 
     &       ampi**4*kp4*p1p2 + ampi**2*kp4*p1p4*p2p3 + 
     &       ampi**2*kp4*p1p3*p2p4 + ame**2*kp3**2*p3p4 + 
     &       kp3**2*p1p2*p3p4 + ame**2*kp3*p3p4**2 + 
     &       kp3*p1p2*p3p4**2 - ame**2*ampi**2*kp3*kp4*1D0 - 
     &       ampi**2*kp3*kp4*p1p2*1D0 - 
     &       ampi**2*kp2*kp4*p1p3*1D0 - 
     &       ampi**2*kp1*kp4*p2p3*1D0 - 
     &       ame**2*ampi**2*kp3*p3p4*1D0 - 
     &       ame**2*ampi**2*kp4*p3p4*1D0 - 
     &       ampi**2*kp3*p1p2*p3p4*1D0 - 
     &       ampi**2*kp4*p1p2*p3p4*1D0 - 
     &       kp2*kp3*p1p3*p3p4*1D0 - kp1*kp3*p2p3*p3p4*1D0 - 
     &       kp3*p1p4*p2p3*p3p4*1D0 - 
     &       kp3*p1p3*p2p4*p3p4*1D0 - 
     &       ampi**2*kp4*p1p3*p2p3*2D0 + 
     &       kp3*p1p3*p2p3*p3p4*2D0 + kp3*kp4*p1p3*p2p3*4D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C22ax45px12pm1p*
     &     (ame**2*ampi**4*kp4 + ampi**4*kp4*p1p2 + 
     &       ampi**2*kp4*p1p4*p2p3 + ampi**2*kp4*p1p3*p2p4 + 
     &       ame**2*kp3**2*p3p4 + kp3**2*p1p2*p3p4 + 
     &       ame**2*kp3*p3p4**2 + kp3*p1p2*p3p4**2 - 
     &       ame**2*ampi**2*kp3*kp4*1D0 - 
     &       ampi**2*kp3*kp4*p1p2*1D0 - 
     &       ampi**2*kp2*kp4*p1p3*1D0 - 
     &       ampi**2*kp1*kp4*p2p3*1D0 - 
     &       ame**2*ampi**2*kp3*p3p4*1D0 - 
     &       ame**2*ampi**2*kp4*p3p4*1D0 - 
     &       ampi**2*kp3*p1p2*p3p4*1D0 - 
     &       ampi**2*kp4*p1p2*p3p4*1D0 - 
     &       kp2*kp3*p1p3*p3p4*1D0 - kp1*kp3*p2p3*p3p4*1D0 - 
     &       kp3*p1p4*p2p3*p3p4*1D0 - 
     &       kp3*p1p3*p2p4*p3p4*1D0 - 
     &       ampi**2*kp4*p1p3*p2p3*2D0 + 
     &       kp3*p1p3*p2p3*p3p4*2D0 + kp3*kp4*p1p3*p2p3*4D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C12ax35x12pm1pp*
     &     (ame**2*ampi**4*kp3 + ame**2*ampi**4*kp4 + 
     &       ampi**4*kp3*p1p2 + ampi**4*kp4*p1p2 + 
     &       ampi**2*kp2*kp3*p1p4 + ampi**2*kp2*kp4*p1p4 + 
     &       ampi**2*kp3*p1p4*p2p3 + ampi**2*kp4*p1p4*p2p3 + 
     &       ampi**2*kp1*kp3*p2p4 + ampi**2*kp1*kp4*p2p4 + 
     &       ampi**2*kp3*p1p3*p2p4 + ampi**2*kp4*p1p3*p2p4 + 
     &       ame**2*kp3*kp4*p3p4 + ame**2*kp4**2*p3p4 + 
     &       kp3*kp4*p1p2*p3p4 + kp4**2*p1p2*p3p4 + 
     &       kp2*kp3*p1p4*p3p4 + kp1*kp3*p2p4*p3p4 + 
     &       ame**2*kp3*p3p4**2 + ame**2*kp4*p3p4**2 + 
     &       kp3*p1p2*p3p4**2 + kp4*p1p2*p3p4**2 - 
     &       ame**2*ampi**2*kp3*kp4*1D0 - 
     &       ame**2*ampi**2*kp4**2*1D0 - 
     &       ampi**2*kp3*kp4*p1p2*1D0 - 
     &       ampi**2*kp4**2*p1p2*1D0 - 
     &       kp2*kp4*p1p4*p3p4*1D0 - kp3*p1p4*p2p3*p3p4*1D0 - 
     &       kp4*p1p4*p2p3*p3p4*1D0 - kp1*kp4*p2p4*p3p4*1D0 - 
     &       kp3*p1p3*p2p4*p3p4*1D0 - 
     &       kp4*p1p3*p2p4*p3p4*1D0 + 
     &       ampi**2*kp1*kp2*kp3*2D0 + 
     &       ame**2*kp3*kp4**2*2D0 + kp3*kp4**2*p1p2*2D0 - 
     &       kp2*kp3*kp4*p1p4*2D0 - 
     &       ampi**2*kp3*p1p3*p2p3*2D0 - 
     &       kp1*kp3*kp4*p2p4*2D0 - 
     &       ampi**2*kp4*p1p4*p2p4*2D0 - 
     &       ame**2*ampi**2*kp3*p3p4*2D0 - 
     &       ame**2*ampi**2*kp4*p3p4*2D0 + 
     &       kp1*kp2*kp4*p3p4*2D0 - 
     &       ampi**2*kp3*p1p2*p3p4*2D0 - 
     &       ampi**2*kp4*p1p2*p3p4*2D0 + 
     &       kp2*kp4*p1p3*p3p4*2D0 + kp1*kp4*p2p3*p3p4*2D0 + 
     &       kp4*p1p3*p2p3*p3p4*2D0 + 
     &       kp3*p1p4*p2p4*p3p4*2D0 - kp3*kp4*p1p4*p2p3*4D0 - 
     &       kp3*kp4*p1p3*p2p4*4D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C12ax45x12pm1pp*
     &     (ame**2*ampi**4*kp3 + ame**2*ampi**4*kp4 + 
     &       ampi**4*kp3*p1p2 + ampi**4*kp4*p1p2 + 
     &       ampi**2*kp2*kp3*p1p3 + ampi**2*kp2*kp4*p1p3 + 
     &       ampi**2*kp1*kp3*p2p3 + ampi**2*kp1*kp4*p2p3 + 
     &       ampi**2*kp3*p1p4*p2p3 + ampi**2*kp4*p1p4*p2p3 + 
     &       ampi**2*kp3*p1p3*p2p4 + ampi**2*kp4*p1p3*p2p4 + 
     &       ame**2*kp3**2*p3p4 + ame**2*kp3*kp4*p3p4 + 
     &       kp3**2*p1p2*p3p4 + kp3*kp4*p1p2*p3p4 + 
     &       kp2*kp4*p1p3*p3p4 + kp1*kp4*p2p3*p3p4 + 
     &       ame**2*kp3*p3p4**2 + ame**2*kp4*p3p4**2 + 
     &       kp3*p1p2*p3p4**2 + kp4*p1p2*p3p4**2 - 
     &       ame**2*ampi**2*kp3**2*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*1D0 - 
     &       ampi**2*kp3**2*p1p2*1D0 - 
     &       ampi**2*kp3*kp4*p1p2*1D0 - 
     &       kp2*kp3*p1p3*p3p4*1D0 - kp1*kp3*p2p3*p3p4*1D0 - 
     &       kp3*p1p4*p2p3*p3p4*1D0 - 
     &       kp4*p1p4*p2p3*p3p4*1D0 - 
     &       kp3*p1p3*p2p4*p3p4*1D0 - 
     &       kp4*p1p3*p2p4*p3p4*1D0 + 
     &       ampi**2*kp1*kp2*kp4*2D0 + 
     &       ame**2*kp3**2*kp4*2D0 + kp3**2*kp4*p1p2*2D0 - 
     &       kp2*kp3*kp4*p1p3*2D0 - kp1*kp3*kp4*p2p3*2D0 - 
     &       ampi**2*kp3*p1p3*p2p3*2D0 - 
     &       ampi**2*kp4*p1p4*p2p4*2D0 - 
     &       ame**2*ampi**2*kp3*p3p4*2D0 + 
     &       kp1*kp2*kp3*p3p4*2D0 - 
     &       ame**2*ampi**2*kp4*p3p4*2D0 - 
     &       ampi**2*kp3*p1p2*p3p4*2D0 - 
     &       ampi**2*kp4*p1p2*p3p4*2D0 + 
     &       kp2*kp3*p1p4*p3p4*2D0 + kp4*p1p3*p2p3*p3p4*2D0 + 
     &       kp1*kp3*p2p4*p3p4*2D0 + kp3*p1p4*p2p4*p3p4*2D0 - 
     &       kp3*kp4*p1p4*p2p3*4D0 - kp3*kp4*p1p3*p2p4*4D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C11apx450pm1p*1D0*
     &     (ame**2*ampi**2*kp3*kp4 + ampi**2*kp3*kp4*p1p2 + 
     &       ampi**2*kp2*kp3*p1p4 + ampi**2*kp1*kp3*p2p4 + 
     &       ame**2*ampi**2*kp3*p3p4 + 
     &       ame**2*ampi**2*kp4*p3p4 + 
     &       ampi**2*kp3*p1p2*p3p4 + ampi**2*kp4*p1p2*p3p4 + 
     &       kp2*kp4*p1p4*p3p4 + kp4*p1p4*p2p3*p3p4 + 
     &       kp1*kp4*p2p4*p3p4 + kp4*p1p3*p2p4*p3p4 - 
     &       ame**2*ampi**4*kp3*1D0 - ampi**4*kp3*p1p2*1D0 - 
     &       ampi**2*kp3*p1p4*p2p3*1D0 - 
     &       ampi**2*kp3*p1p3*p2p4*1D0 - 
     &       ame**2*kp4**2*p3p4*1D0 - kp4**2*p1p2*p3p4*1D0 - 
     &       ame**2*kp4*p3p4**2*1D0 - kp4*p1p2*p3p4**2*1D0 + 
     &       ampi**2*kp3*p1p4*p2p4*2D0 - 
     &       kp4*p1p4*p2p4*p3p4*2D0 - kp3*kp4*p1p4*p2p4*4D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C22ax35px12pm1p*1D0*
     &     (ame**2*ampi**2*kp3*kp4 + ampi**2*kp3*kp4*p1p2 + 
     &       ampi**2*kp2*kp3*p1p4 + ampi**2*kp1*kp3*p2p4 + 
     &       ame**2*ampi**2*kp3*p3p4 + 
     &       ame**2*ampi**2*kp4*p3p4 + 
     &       ampi**2*kp3*p1p2*p3p4 + ampi**2*kp4*p1p2*p3p4 + 
     &       kp2*kp4*p1p4*p3p4 + kp4*p1p4*p2p3*p3p4 + 
     &       kp1*kp4*p2p4*p3p4 + kp4*p1p3*p2p4*p3p4 - 
     &       ame**2*ampi**4*kp3*1D0 - ampi**4*kp3*p1p2*1D0 - 
     &       ampi**2*kp3*p1p4*p2p3*1D0 - 
     &       ampi**2*kp3*p1p3*p2p4*1D0 - 
     &       ame**2*kp4**2*p3p4*1D0 - kp4**2*p1p2*p3p4*1D0 - 
     &       ame**2*kp4*p3p4**2*1D0 - kp4*p1p2*p3p4**2*1D0 + 
     &       ampi**2*kp3*p1p4*p2p4*2D0 - 
     &       kp4*p1p4*p2p4*p3p4*2D0 - kp3*kp4*p1p4*p2p4*4D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C22apx45x12pm1p*1D0*
     &     (ame**2*ampi**2*kp3**2 + ampi**2*kp3**2*p1p2 + 
     &       ame**2*ampi**2*kp3*p3p4 + 
     &       ame**2*ampi**2*kp4*p3p4 + 
     &       ampi**2*kp3*p1p2*p3p4 + ampi**2*kp4*p1p2*p3p4 + 
     &       kp2*kp4*p1p3*p3p4 + kp1*kp4*p2p3*p3p4 + 
     &       kp4*p1p4*p2p3*p3p4 + kp4*p1p3*p2p4*p3p4 - 
     &       ame**2*ampi**4*kp3*1D0 - ampi**4*kp3*p1p2*1D0 - 
     &       ampi**2*kp2*kp3*p1p3*1D0 - 
     &       ampi**2*kp1*kp3*p2p3*1D0 - 
     &       ampi**2*kp3*p1p4*p2p3*1D0 - 
     &       ampi**2*kp3*p1p3*p2p4*1D0 - 
     &       ame**2*kp3*kp4*p3p4*1D0 - 
     &       kp3*kp4*p1p2*p3p4*1D0 - ame**2*kp4*p3p4**2*1D0 - 
     &       kp4*p1p2*p3p4**2*1D0 + ampi**2*kp1*kp2*kp3*2D0 + 
     &       ame**2*kp3*kp4**2*2D0 + kp3*kp4**2*p1p2*2D0 + 
     &       ampi**2*kp2*kp3*p1p4*2D0 - 
     &       kp2*kp3*kp4*p1p4*2D0 + 
     &       ampi**2*kp1*kp3*p2p4*2D0 - 
     &       kp1*kp3*kp4*p2p4*2D0 + 
     &       ampi**2*kp3*p1p4*p2p4*2D0 + 
     &       kp1*kp2*kp4*p3p4*2D0 - kp4*p1p4*p2p4*p3p4*2D0 - 
     &       kp3*kp4*p1p4*p2p4*4D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C11ax12x340ppp*1D0*
     &     (ame**2*ampi**2*kp3**2 + ame**2*ampi**2*kp4**2 + 
     &       ampi**2*kp3**2*p1p2 + ampi**2*kp4**2*p1p2 + 
     &       ame**2*kp3**2*p3p4 + ame**2*kp4**2*p3p4 + 
     &       kp3**2*p1p2*p3p4 + kp4**2*p1p2*p3p4 - 
     &       ampi**2*kp2*kp3*p1p3*1D0 - 
     &       ampi**2*kp2*kp4*p1p3*1D0 - 
     &       ampi**2*kp2*kp3*p1p4*1D0 - 
     &       ampi**2*kp2*kp4*p1p4*1D0 - 
     &       ampi**2*kp1*kp3*p2p3*1D0 - 
     &       ampi**2*kp1*kp4*p2p3*1D0 - 
     &       ampi**2*kp1*kp3*p2p4*1D0 - 
     &       ampi**2*kp1*kp4*p2p4*1D0 - 
     &       kp2*kp3*p1p3*p3p4*1D0 - kp2*kp4*p1p3*p3p4*1D0 - 
     &       kp2*kp3*p1p4*p3p4*1D0 - kp2*kp4*p1p4*p3p4*1D0 - 
     &       kp1*kp3*p2p3*p3p4*1D0 - kp1*kp4*p2p3*p3p4*1D0 - 
     &       kp1*kp3*p2p4*p3p4*1D0 - kp1*kp4*p2p4*p3p4*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*2D0 - 
     &       ampi**2*kp3*kp4*p1p2*2D0 + 
     &       ampi**2*kp3*p1p3*p2p3*2D0 - 
     &       ampi**2*kp4*p1p3*p2p3*2D0 - 
     &       ampi**2*kp3*p1p4*p2p4*2D0 + 
     &       ampi**2*kp4*p1p4*p2p4*2D0 - 
     &       ame**2*kp3*kp4*p3p4*2D0 - 
     &       kp3*kp4*p1p2*p3p4*2D0 + kp3*p1p3*p2p3*p3p4*2D0 - 
     &       kp4*p1p3*p2p3*p3p4*2D0 - 
     &       kp3*p1p4*p2p4*p3p4*2D0 + 
     &       kp4*p1p4*p2p4*p3p4*2D0 + kp3*kp4*p1p3*p2p3*4D0 + 
     &       kp3*kp4*p1p4*p2p3*4D0 + kp3*kp4*p1p3*p2p4*4D0 + 
     &       kp3*kp4*p1p4*p2p4*4D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D0ap0x12px45x34m1ppp*
     &     (-(kp4*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(kp3*(ame**2 + p1p2)*(ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(ame**2*kp4**2 + kp4**2*p1p2 + 
     &          ampi**2*kp3*(ame**2 + p1p2) + kp4*p1p4*p2p3 - 
     &          kp2*kp4*p1p4*1D0 - ampi**2*p1p4*p2p3*1D0 - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0) + 
     &             ampi**2*(kp2 - p2p3*1D0)) + 
     &          p2p4*(p1p3*(kp4 - ampi**2*1D0) + 
     &             p1p4*(ampi**2 - kp4*2D0))))*
     &     (m1**2 + p3p4*4D0))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D22ap0x12px35x34m1ppp*
     &     (-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(p1p4*(kp2 + p2p3) + p1p3*p2p4 + 
     &             (ame**2 + p1p2)*(ampi**2 - kp4*1D0) - 
     &             p1p3*p2p3*2D0 + kp1*(p2p4 + kp2*2D0))) + 
     &       kp4*(kp3**2*x12 - 
     &          kp3*p1p3*2D0*(kp2 + p2p3*2D0) + 
     &          kp1*(-(ampi**2*p2p4*1D0) + ampi**2*kp2*2D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)*2D0) - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + p1p3*p2p4 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             (kp2 + p2p3)*1D0*(-(p1p4*1D0) + p1p3*2D0)))
     &       )*(m1**2 + p3p4*4D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D0ap0x12px35x34m1ppp*
     &     (-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(kp3**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*kp2 - ampi**2*p2p4*1D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)) + 
     &          ampi**2*(kp4*(ame**2 + p1p2) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &             (p1p3 - p1p4*1D0)*(kp2 + p2p3 - p2p4*1D0))
     &           + kp3*(p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             p1p3*(p2p4 - p2p3*2D0))))*
     &     (m1**2 + p3p4*4D0))/(kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (D23ap0x12px45x34m1ppp*1D0*
     &     (-(kp3**2*(ame**2 + p1p2)*1D0*
     &          (ampi**2 - p3p4*1D0 - kp4*2D0)) + 
     &       kp4*((ame**2 + p1p2)*p3p4**2 + 
     &          p3p4*(kp2*p1p3 - p1p3*p2p4*1D0 - 
     &             ame**2*ampi**2*2D0 - ampi**2*p1p2*2D0 + 
     &             p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &          ampi**2*(ampi**2*(ame**2 + p1p2) + kp2*p1p3 + 
     &             p1p4*p2p3 + kp1*(p2p3 + kp2*2D0) + 
     &             p2p4*(p1p3 - p1p4*2D0))) + 
     &       kp3*((ame**2 + p1p2)*p3p4**2 + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             kp2*p1p3 + p1p3*p2p4 + 
     &             p2p3*(kp1 + p1p4 - p1p3*2D0)) - 
     &          p3p4*1D0*
     &           (p1p3*p2p4 + ampi**2*x12 + 
     &             kp2*(p1p3 - p1p4*2D0) + 
     &             p1p4*(p2p3 - p2p4*2D0) + 
     &             kp1*(p2p3 - kp2*2D0 - p2p4*2D0)) + 
     &          kp4*(ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) - 
     &             p2p3*2D0*(kp1 + p1p4*2D0) - 
     &             p1p3*2D0*(kp2 + p2p4*2D0))))*
     &     (m1**2 + p3p4*4D0))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D22ap0x12px45x34m1ppp*
     &     (ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*p3p4*(kp2*p1p3 + p1p4*p2p3 + 
     &          (ame**2 + p1p2)*(ampi**2 - p3p4*1D0) + 
     &          kp1*(p2p3 + kp2*2D0) + p2p4*(p1p3 - p1p4*2D0))
     &         + kp3*(ampi**2*(ame**2 + p1p2)*p3p4 + 
     &          kp4**2*x12 + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0)*2D0 + 
     &             ampi**2*(-(p2p3*1D0) + kp2*2D0)) - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + p1p3*p2p4 + 
     &             kp2*(p1p3 - p1p4*2D0) + 
     &             p1p4*(p2p3 - p2p4*2D0)) + 
     &          kp4*(-((ame**2 + p1p2)*p3p4*1D0) - 
     &             p1p4*2D0*(kp2 + p2p4*2D0))))*
     &     (m1**2 + p3p4*4D0))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D00ap0x12px35x34m1ppp*
     &     (kp4*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0)))*
     &     (m1**2 + p3p4*4D0))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D00ap0x12px45x34m1ppp*
     &     (kp4*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0)))*
     &     (m1**2 + p3p4*4D0))/(kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C0appx34pm1p*5.D-1*
     &     (kp4*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0)))*
     &     (m1**2 + p3p4*4D0))/(kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (D3ap0x12px45x34m1ppp*1D0*
     &     (-(kp3**2*(ame**2 + p1p2)*1D0*
     &          (ampi**2 - kp4*1D0 - p3p4*1D0)) + 
     &       kp4*((ame**2 + p1p2)*p3p4**2 + 
     &          ampi**2*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)) + 
     &          p3p4*(kp2*p1p3 - p1p3*p2p4*1D0 - 
     &             ame**2*ampi**2*2D0 - ampi**2*p1p2*2D0 + 
     &             p2p3*(kp1 - p1p4*1D0 + p1p3*2D0))) + 
     &       kp3*(kp4*(ame**2 + p1p2)*p3p4 + 
     &          (ame**2 + p1p2)*p3p4**2 + 
     &          p3p4*((kp1 + p1p4 - p1p3*1D0)*
     &              (kp2 + p2p4 - p2p3*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             kp2*p1p3 + p1p3*p2p4 + 
     &             p2p3*(kp1 + p1p4 - p1p3*2D0)) - 
     &          kp4*1D0*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             kp2*p1p3 + p1p3*p2p4*3D0 + 
     &             p2p3*(kp1 - p1p3*2D0 + p1p4*3D0))))*
     &     (m1**2 + p3p4*4D0))/(kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (D3ap0x12px35x34m1ppp*1D0*
     &     (kp4*(ame**2*p3p4**2 - 
     &          kp4*(ame**2 + p1p2)*1D0*
     &           (ampi**2 - p3p4*1D0) + 
     &          p1p2*(ampi**2 - p3p4*1D0)**2 + 
     &          p3p4*((kp1 + p1p3 - p1p4*1D0)*
     &              (kp2 + p2p3 - p2p4*1D0) - 
     &             ame**2*ampi**2*2D0) + 
     &          ampi**2*(ame**2*ampi**2 + (kp1 + p1p3)*p2p4 + 
     &             p1p4*(kp2 + p2p3 - p2p4*2D0))) + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          kp4**2*(ame**2 + p1p2) + 
     &          (ame**2 + p1p2)*p3p4**2 - 
     &          ampi**2*1D0*(p1p3 - p1p4*1D0)*
     &           (p2p3 - p2p4*1D0) + 
     &          kp1*(ampi**2*kp2 + p2p4*(p3p4 - kp4*1D0)) - 
     &          p3p4*1D0*
     &           (p1p3*p2p4 + ame**2*ampi**2*2D0 + 
     &             ampi**2*p1p2*2D0 + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)) + 
     &          kp4*(ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) - 
     &             p1p3*p2p4*3D0 - 
     &             p1p4*1D0*(kp2 - p2p4*2D0 + p2p3*3D0))))*
     &     (m1**2 + p3p4*4D0))/(kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (D13ap0x12px35x34m1ppp*1D0*(m1**2 + p3p4*4D0)*
     &     (-(ampi**2*kp3**2*(ame**2 + p1p2)*1D0) + 
     &       kp4*((ame**2 + p1p2)*p3p4**2 + 
     &          p3p4*(kp2*p1p3 - p1p3*p2p4*1D0 - 
     &             ame**2*ampi**2*2D0 - ampi**2*p1p2*2D0 + 
     &             p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             (kp1 + p1p3)*p2p4 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p1p4*(kp2 + p2p3 - p2p4*2D0))) + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 + ampi**2*p1p3*p2p4 + 
     &          (ame**2 + p1p2)*p3p4**2 - 
     &          ampi**2*p1p3*p2p3*2D0 + 
     &          p3p4*(kp4*x12 + p2p4*(kp1 - p1p3*1D0) - 
     &             ame**2*ampi**2*2D0 - ampi**2*p1p2*2D0 + 
     &             p1p4*(kp2 - p2p3*1D0 + p2p4*2D0)) - 
     &          kp4*p1p4*p2p3*4D0 - kp4*p1p3*p2p4*4D0)))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (D13ap0x12px45x34m1ppp*1D0*(m1**2 + p3p4*4D0)*
     &     (-(ampi**2*kp3**2*(ame**2 + p1p2)*1D0) + 
     &       kp4*((ame**2 + p1p2)*p3p4**2 + 
     &          p3p4*(kp2*p1p3 - p1p3*p2p4*1D0 - 
     &             ame**2*ampi**2*2D0 - ampi**2*p1p2*2D0 + 
     &             p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             (kp1 + p1p3)*p2p4 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p1p4*(kp2 + p2p3 - p2p4*2D0))) + 
     &       kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &          ampi**2*kp2*p1p3 + ampi**2*kp1*p2p3 + 
     &          ampi**2*p1p4*p2p3 + ampi**2*p1p3*p2p4 + 
     &          (ame**2 + p1p2)*p3p4**2 - 
     &          ampi**2*p1p3*p2p3*2D0 + 
     &          p3p4*(kp4*x12 + p2p4*(kp1 - p1p3*1D0) - 
     &             ame**2*ampi**2*2D0 - ampi**2*p1p2*2D0 + 
     &             p1p4*(kp2 - p2p3*1D0 + p2p4*2D0)) - 
     &          kp4*p1p4*p2p3*4D0 - kp4*p1p3*p2p4*4D0)))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C11apx12x45m1pp*
     &     (kp3*p3p4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &          kp2*p1p3 + p1p3*p2p4 - 
     &          (ame**2 + p1p2)*(kp3 + p3p4)*1D0 + 
     &          p2p3*(kp1 + p1p4 - p1p3*2D0)) + 
     &       kp4*(ampi**2*
     &           (kp2*p1p3 - p1p3*p2p4*1D0 + 
     &             ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &          kp3*(ampi**2*(ame**2 + p1p2) - p1p3*p2p3*4D0))
     &       ))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C11apx350pm1p*(kp3*p3p4*
     &        (ame**2*ampi**2 + ampi**2*p1p2 + kp2*p1p3 + 
     &          p1p3*p2p4 - 
     &          (ame**2 + p1p2)*(kp3 + p3p4)*1D0 + 
     &          p2p3*(kp1 + p1p4 - p1p3*2D0)) + 
     &       kp4*(ampi**2*
     &           (kp2*p1p3 - p1p3*p2p4*1D0 + 
     &             ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &          kp3*(ampi**2*(ame**2 + p1p2) - p1p3*p2p3*4D0))
     &       ))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C1apx350pm1p*(kp3*p3p4*
     &        ((ame**2 + p1p2)*(kp3 + p3p4) - 
     &          ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 - 
     &          kp2*p1p3*1D0 - p1p3*p2p4*1D0 - 
     &          p2p3*1D0*(kp1 + p1p4 - p1p3*2D0)) + 
     &       kp4*(ampi**2*
     &           (p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)) - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*4D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C0apx12x45m1pp*2D0*
     &     (kp3*p3p4*((ame**2 + p1p2)*(kp3 + p3p4) - 
     &          ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 - 
     &          kp2*p1p3*1D0 - p1p3*p2p4*1D0 - 
     &          p2p3*1D0*(kp1 + p1p4 - p1p3*2D0)) + 
     &       kp4*(ampi**2*
     &           (p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)) - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*4D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C1apx12x45m1pp*3D0*
     &     (kp3*p3p4*((ame**2 + p1p2)*(kp3 + p3p4) - 
     &          ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 - 
     &          kp2*p1p3*1D0 - p1p3*p2p4*1D0 - 
     &          p2p3*1D0*(kp1 + p1p4 - p1p3*2D0)) + 
     &       kp4*(ampi**2*
     &           (p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)) - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*4D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (D11ap0x12px35x34m1ppp*1D0*(m1**2 + p3p4*4D0)*
     &     (kp3*p3p4*((ame**2 + p1p2)*(kp3 + p3p4) - 
     &          ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 - 
     &          kp2*p1p3*1D0 - p1p3*p2p4*1D0 - 
     &          p2p3*1D0*(kp1 + p1p4 - p1p3*2D0)) + 
     &       kp4*(ampi**2*
     &           (p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)) - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*4D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (D33ap0x12px45x34m1ppp*1D0*(m1**2 + p3p4*4D0)*
     &     (kp3*p3p4*((ame**2 + p1p2)*(kp3 + p3p4) - 
     &          ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 - 
     &          kp2*p1p3*1D0 - p1p3*p2p4*1D0 - 
     &          p2p3*1D0*(kp1 + p1p4 - p1p3*2D0)) + 
     &       kp4*(ampi**2*
     &           (p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)) - 
     &          kp3*1D0*(ampi**2*(ame**2 + p1p2) - 
     &             p1p3*p2p3*4D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (D23ap0x12px35x34m1ppp*1D0*(m1**2 + p3p4*4D0)*
     &     (kp4*((ame**2 + p1p2)*p3p4**2 - 
     &          p3p4*1D0*
     &           ((kp1 + p1p3)*p2p4 + 
     &             (ame**2 + p1p2)*
     &              (-(kp4*1D0) + ampi**2*2D0) - 
     &             (kp2 + p2p3)*1D0*
     &              (-(p1p4*1D0) + (kp1 + p1p3)*2D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             (kp1 + p1p3)*p2p4 - 
     &             kp4*(ame**2 + p1p2)*1D0 + 
     &             p1p4*(kp2 + p2p3 - p2p4*2D0))) + 
     &       kp3*((ame**2 + p1p2)*p3p4**2 + kp4**2*x12 + 
     &          kp1*(ampi**2*kp2*2D0 + 
     &             p2p4*(ampi**2 + p3p4 - kp4*2D0)) + 
     &          ampi**2*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &             p1p4*(kp2 + p2p3) + p1p3*(p2p4 - p2p3*2D0))
     &            - p3p4*1D0*
     &           (p1p3*p2p4 + ame**2*ampi**2*2D0 + 
     &             ampi**2*p1p2*2D0 + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)) + 
     &          kp4*(ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) - 
     &             p1p4*2D0*(kp2 + p2p3*2D0) - p1p3*p2p4*4D0))
     &       ))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C11apx12x35m1pp*
     &     (kp4*p3p4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &          (kp1 + p1p3)*p2p4 - 
     &          (ame**2 + p1p2)*(kp4 + p3p4)*1D0 + 
     &          p1p4*(kp2 + p2p3 - p2p4*2D0)) + 
     &       kp3*(ampi**2*
     &           (ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             p2p4*(kp1 - p1p3*1D0) + 
     &             p1p4*(kp2 - p2p3*1D0 + p2p4*2D0)) + 
     &          kp4*(ampi**2*(ame**2 + p1p2) - p1p4*p2p4*4D0))
     &       ))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C11apx450pm1p*(kp4*p3p4*
     &        (ame**2*ampi**2 + ampi**2*p1p2 + 
     &          (kp1 + p1p3)*p2p4 - 
     &          (ame**2 + p1p2)*(kp4 + p3p4)*1D0 + 
     &          p1p4*(kp2 + p2p3 - p2p4*2D0)) + 
     &       kp3*(ampi**2*
     &           (ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             p2p4*(kp1 - p1p3*1D0) + 
     &             p1p4*(kp2 - p2p3*1D0 + p2p4*2D0)) + 
     &          kp4*(ampi**2*(ame**2 + p1p2) - p1p4*p2p4*4D0))
     &       ))/(kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C1apx450pm1p*1D0*
     &     (kp4*p3p4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &          (kp1 + p1p3)*p2p4 - 
     &          (ame**2 + p1p2)*(kp4 + p3p4)*1D0 + 
     &          p1p4*(kp2 + p2p3 - p2p4*2D0)) + 
     &       kp3*(ampi**2*
     &           (ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             p2p4*(kp1 - p1p3*1D0) + 
     &             p1p4*(kp2 - p2p3*1D0 + p2p4*2D0)) + 
     &          kp4*(ampi**2*(ame**2 + p1p2) - p1p4*p2p4*4D0))
     &       ))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C0apx12x35m1pp*2D0*
     &     (kp4*p3p4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &          (kp1 + p1p3)*p2p4 - 
     &          (ame**2 + p1p2)*(kp4 + p3p4)*1D0 + 
     &          p1p4*(kp2 + p2p3 - p2p4*2D0)) + 
     &       kp3*(ampi**2*
     &           (ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             p2p4*(kp1 - p1p3*1D0) + 
     &             p1p4*(kp2 - p2p3*1D0 + p2p4*2D0)) + 
     &          kp4*(ampi**2*(ame**2 + p1p2) - p1p4*p2p4*4D0))
     &       ))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C1apx12x35m1pp*3D0*
     &     (kp4*p3p4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &          (kp1 + p1p3)*p2p4 - 
     &          (ame**2 + p1p2)*(kp4 + p3p4)*1D0 + 
     &          p1p4*(kp2 + p2p3 - p2p4*2D0)) + 
     &       kp3*(ampi**2*
     &           (ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             p2p4*(kp1 - p1p3*1D0) + 
     &             p1p4*(kp2 - p2p3*1D0 + p2p4*2D0)) + 
     &          kp4*(ampi**2*(ame**2 + p1p2) - p1p4*p2p4*4D0))
     &       ))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D11ap0x12px45x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (kp4*p3p4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &          (kp1 + p1p3)*p2p4 - 
     &          (ame**2 + p1p2)*(kp4 + p3p4)*1D0 + 
     &          p1p4*(kp2 + p2p3 - p2p4*2D0)) + 
     &       kp3*(ampi**2*
     &           (ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             p2p4*(kp1 - p1p3*1D0) + 
     &             p1p4*(kp2 - p2p3*1D0 + p2p4*2D0)) + 
     &          kp4*(ampi**2*(ame**2 + p1p2) - p1p4*p2p4*4D0))
     &       ))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D33ap0x12px35x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (kp4*p3p4*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &          (kp1 + p1p3)*p2p4 - 
     &          (ame**2 + p1p2)*(kp4 + p3p4)*1D0 + 
     &          p1p4*(kp2 + p2p3 - p2p4*2D0)) + 
     &       kp3*(ampi**2*
     &           (ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) + 
     &             p2p4*(kp1 - p1p3*1D0) + 
     &             p1p4*(kp2 - p2p3*1D0 + p2p4*2D0)) + 
     &          kp4*(ampi**2*(ame**2 + p1p2) - p1p4*p2p4*4D0))
     &       ))/(kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C11ax45x12pm1pp*
     &     (-(kp3*1D0*(p1p4 + p2p4 - p3p4*1D0)*
     &          (ame**2*ampi**2 + ampi**2*p1p2 + kp2*p1p3 + 
     &            p1p3*p2p4 - 
     &            (ame**2 + p1p2)*(kp3 + p3p4)*1D0 + 
     &            p2p3*(kp1 + p1p4 - p1p3*2D0))) + 
     &       kp4*(-(1D0*(ampi**2 - p1p3*1D0 - p2p3*1D0)*
     &             (p1p3*p2p4 - kp2*p1p3*1D0 + 
     &               ame**2*(ampi**2 - p3p4*1D0) + 
     &               p1p2*(ampi**2 - p3p4*1D0) + 
     &               p2p3*(p1p4 - kp1*1D0 - p1p3*2D0))) + 
     &          kp3*((ame**2 + p1p2)*(ampi**2 + p1p3) + 
     &             p2p3*(ame**2 + p1p2 - p1p3*4D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C11ax35x12pm1pp*
     &     (-(kp4*1D0*(p1p3 + p2p3 - p3p4*1D0)*
     &          (ame**2*ampi**2 + ampi**2*p1p2 + 
     &            (kp1 + p1p3)*p2p4 - 
     &            (ame**2 + p1p2)*(kp4 + p3p4)*1D0 + 
     &            p1p4*(kp2 + p2p3 - p2p4*2D0))) + 
     &       kp3*(-(1D0*(ampi**2 - p1p4*1D0 - p2p4*1D0)*
     &             (p2p4*(p1p3 - kp1*1D0) + 
     &               ame**2*(ampi**2 - p3p4*1D0) + 
     &               p1p2*(ampi**2 - p3p4*1D0) + 
     &               p1p4*(p2p3 - kp2*1D0 - p2p4*2D0))) + 
     &          kp4*((ame**2 + p1p2)*(ampi**2 + p1p4) + 
     &             p2p4*(ame**2 + p1p2 - p1p4*4D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D12ap0x12px35x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (-(kp3*(ame**2 + p1p2)*p3p4**2*2D0) + 
     &       p3p4*(ampi**2*kp4*x12 - 
     &          kp3**2*(ame**2 + p1p2)*1D0 + 
     &          kp3*(kp2*(p1p3 + p1p4) + p1p3*p2p4*2D0 + 
     &             (ame**2 + p1p2)*
     &              (-(kp4*1D0) + ampi**2*2D0) + 
     &             kp1*(p2p3 + p2p4 + kp2*2D0) + 
     &             p2p3*2D0*(p1p4 - p1p3*2D0))) + 
     &       kp4*(kp3**2*x12 + 
     &          kp1*(-(ampi**2*p2p4*1D0) + ampi**2*kp2*2D0 + 
     &             p2p3*(-(kp3*2D0) + ampi**2*3D0)) + 
     &          ampi**2*(kp4*(ame**2 + p1p2) - 
     &             ampi**2*p1p2*2D0 - 
     &             2D0*(ame**2*ampi**2 + p1p4*p2p3 + 
     &                p1p3*(p2p4 - p2p3*2D0)) + 
     &             kp2*(-(p1p4*1D0) + p1p3*3D0)) + 
     &          kp3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*2D0*(kp2 + p2p3*4D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D12ap0x12px45x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*p3p4*(kp2*(p1p3 + p1p4) + p1p4*p2p3*2D0 + 
     &          kp1*(p2p3 + p2p4 + kp2*2D0) + 
     &          p2p4*2D0*(p1p3 - p1p4*2D0) + 
     &          (ame**2 + p1p2)*
     &           (-(kp4*1D0) + ampi**2*2D0 - p3p4*2D0)) + 
     &       kp3*(kp4**2*x12 + 
     &          kp1*(ampi**2*(-(p2p3*1D0) + kp2*2D0) + 
     &             p2p4*(-(kp4*2D0) + ampi**2*3D0)) + 
     &          ampi**2*(p1p2*(p3p4 - ampi**2*1D0)*2D0 - 
     &             2D0*(p1p3*p2p4 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p4*(p2p3 - p2p4*2D0)) - 
     &             kp2*1D0*(p1p3 - p1p4*3D0)) + 
     &          kp4*(ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) - 
     &             p1p4*2D0*(kp2 + p2p4*4D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (B0ax35pm1*2.5D-1*
     &     (-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(p1p4*(kp2 + p2p3) + p1p3*p2p4 + 
     &             (ame**2 + p1p2)*(ampi**2 - kp4*1D0) - 
     &             p1p3*p2p3*2D0 + kp1*(p2p4 + kp2*2D0))) + 
     &       kp4*(kp3**2*x12 - 
     &          kp3*p1p3*2D0*(kp2 + p2p3*2D0) + 
     &          kp1*(-(ampi**2*p2p4*1D0) + ampi**2*kp2*2D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)*2D0) - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + p1p3*p2p4 - 
     &             kp4*(ame**2 + p1p2)*1D0 - 
     &             (kp2 + p2p3)*1D0*(-(p1p4*1D0) + p1p3*2D0)))
     &       )*(-(m1**2*1D0) + ampi**2*4D0 + kp3*6D0))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2*x35) + 
     &  (B0ax35pm1*2.5D-1*
     &     (-(kp3*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          kp3*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - kp4*(ame**2 + p1p2)*1D0 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp4*(kp3**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*kp2 - ampi**2*p2p4*1D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)) + 
     &          ampi**2*(kp4*(ame**2 + p1p2) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &             (p1p3 - p1p4*1D0)*(kp2 + p2p3 - p2p4*1D0))
     &           + kp3*(p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             p1p3*(p2p4 - p2p3*2D0))))*
     &     (-(m1**2*1D0) + ampi**2*4D0 + kp3*6D0))/
     &   (kp3**2*kp4*(ame**2 + p1p2)**2*x35) + 
     &  (B0ax45pm1*2.5D-1*
     &     (-(kp4*(ame**2 + p1p2)*p3p4**2*1D0) + 
     &       p3p4*(kp3*(ame**2 + p1p2)*(ampi**2 - kp4*1D0) + 
     &          kp4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)))
     &        + kp3*(ame**2*kp4**2 + kp4**2*p1p2 + 
     &          ampi**2*kp3*(ame**2 + p1p2) + kp4*p1p4*p2p3 - 
     &          kp2*kp4*p1p4*1D0 - ampi**2*p1p4*p2p3*1D0 - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 - 
     &             p1p3*p2p3*1D0 + kp2*(p1p3 - p1p4*1D0)) + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0) + 
     &             ampi**2*(kp2 - p2p3*1D0)) + 
     &          p2p4*(p1p3*(kp4 - ampi**2*1D0) + 
     &             p1p4*(ampi**2 - kp4*2D0))))*
     &     (-(m1**2*1D0) + ampi**2*4D0 + kp4*6D0))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2*x45) + 
     &  (B0ax45pm1*2.5D-1*
     &     (ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*p3p4*(kp2*p1p3 + p1p4*p2p3 + 
     &          (ame**2 + p1p2)*(ampi**2 - p3p4*1D0) + 
     &          kp1*(p2p3 + kp2*2D0) + p2p4*(p1p3 - p1p4*2D0))
     &         + kp3*(ampi**2*(ame**2 + p1p2)*p3p4 + 
     &          kp4**2*x12 + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0)*2D0 + 
     &             ampi**2*(-(p2p3*1D0) + kp2*2D0)) - 
     &          ampi**2*1D0*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + p1p3*p2p4 + 
     &             kp2*(p1p3 - p1p4*2D0) + 
     &             p1p4*(p2p3 - p2p4*2D0)) + 
     &          kp4*(-((ame**2 + p1p2)*p3p4*1D0) - 
     &             p1p4*2D0*(kp2 + p2p4*2D0))))*
     &     (-(m1**2*1D0) + ampi**2*4D0 + kp4*6D0))/
     &   (kp3*kp4**2*(ame**2 + p1p2)**2*x45) - 
     &  (C00apx350pm1p*1D0*
     &     (kp2*kp4*p1p3 + kp2*kp3*p1p4 + kp1*kp4*p2p3 + 
     &       kp1*kp3*p2p4 + ame**2*kp3*p3p4 + 
     &       ame**2*kp4*p3p4 + kp3*p1p2*p3p4 + 
     &       kp4*p1p2*p3p4 - ame**2*ampi**2*kp3*1D0 - 
     &       ame**2*ampi**2*kp4*1D0 - ampi**2*kp3*p1p2*1D0 - 
     &       ampi**2*kp4*p1p2*1D0 - kp3*p1p4*p2p3*1D0 - 
     &       kp4*p1p4*p2p3*1D0 - kp3*p1p3*p2p4*1D0 - 
     &       kp4*p1p3*p2p4*1D0 + kp3*kp4*p1p2*2D0 + 
     &       kp4*p1p3*p2p3*2D0 + kp3*p1p4*p2p4*2D0 + 
     &       ame**2*kp3*kp4*6D0))/(kp3*kp4*(ame**2 + p1p2)**2)
     &    - (C00apx35x12pm1p*1D0*
     &     (kp2*kp4*p1p3 + kp2*kp3*p1p4 + kp1*kp4*p2p3 + 
     &       kp1*kp3*p2p4 + ame**2*kp3*p3p4 + 
     &       ame**2*kp4*p3p4 + kp3*p1p2*p3p4 + 
     &       kp4*p1p2*p3p4 - ame**2*ampi**2*kp3*1D0 - 
     &       ame**2*ampi**2*kp4*1D0 - ampi**2*kp3*p1p2*1D0 - 
     &       ampi**2*kp4*p1p2*1D0 - kp3*p1p4*p2p3*1D0 - 
     &       kp4*p1p4*p2p3*1D0 - kp3*p1p3*p2p4*1D0 - 
     &       kp4*p1p3*p2p4*1D0 + kp3*kp4*p1p2*2D0 + 
     &       kp4*p1p3*p2p3*2D0 + kp3*p1p4*p2p4*2D0 + 
     &       ame**2*kp3*kp4*6D0))/(kp3*kp4*(ame**2 + p1p2)**2)
     &    - (C00apx450pm1p*1D0*
     &     (kp2*kp4*p1p3 + kp2*kp3*p1p4 + kp1*kp4*p2p3 + 
     &       kp1*kp3*p2p4 + ame**2*kp3*p3p4 + 
     &       ame**2*kp4*p3p4 + kp3*p1p2*p3p4 + 
     &       kp4*p1p2*p3p4 - ame**2*ampi**2*kp3*1D0 - 
     &       ame**2*ampi**2*kp4*1D0 - ampi**2*kp3*p1p2*1D0 - 
     &       ampi**2*kp4*p1p2*1D0 - kp3*p1p4*p2p3*1D0 - 
     &       kp4*p1p4*p2p3*1D0 - kp3*p1p3*p2p4*1D0 - 
     &       kp4*p1p3*p2p4*1D0 + kp3*kp4*p1p2*2D0 + 
     &       kp4*p1p3*p2p3*2D0 + kp3*p1p4*p2p4*2D0 + 
     &       ame**2*kp3*kp4*6D0))/(kp3*kp4*(ame**2 + p1p2)**2)
     &    - (C00apx45x12pm1p*1D0*
     &     (kp2*kp4*p1p3 + kp2*kp3*p1p4 + kp1*kp4*p2p3 + 
     &       kp1*kp3*p2p4 + ame**2*kp3*p3p4 + 
     &       ame**2*kp4*p3p4 + kp3*p1p2*p3p4 + 
     &       kp4*p1p2*p3p4 - ame**2*ampi**2*kp3*1D0 - 
     &       ame**2*ampi**2*kp4*1D0 - ampi**2*kp3*p1p2*1D0 - 
     &       ampi**2*kp4*p1p2*1D0 - kp3*p1p4*p2p3*1D0 - 
     &       kp4*p1p4*p2p3*1D0 - kp3*p1p3*p2p4*1D0 - 
     &       kp4*p1p3*p2p4*1D0 + kp3*kp4*p1p2*2D0 + 
     &       kp4*p1p3*p2p3*2D0 + kp3*p1p4*p2p4*2D0 + 
     &       ame**2*kp3*kp4*6D0))/(kp3*kp4*(ame**2 + p1p2)**2)
     &    + (C00a0x12x34ppp*2D0*
     &     (kp2*kp4*p1p3 + kp2*kp3*p1p4 + kp1*kp4*p2p3 + 
     &       kp1*kp3*p2p4 + ame**2*kp3*p3p4 + 
     &       ame**2*kp4*p3p4 + kp3*p1p2*p3p4 + 
     &       kp4*p1p2*p3p4 - ame**2*ampi**2*kp3*1D0 - 
     &       ame**2*ampi**2*kp4*1D0 - ampi**2*kp3*p1p2*1D0 - 
     &       ampi**2*kp4*p1p2*1D0 - kp3*p1p4*p2p3*1D0 - 
     &       kp4*p1p4*p2p3*1D0 - kp3*p1p3*p2p4*1D0 - 
     &       kp4*p1p3*p2p4*1D0 + kp3*kp4*p1p2*2D0 + 
     &       kp4*p1p3*p2p3*2D0 + kp3*p1p4*p2p4*2D0 + 
     &       ame**2*kp3*kp4*6D0))/(kp3*kp4*(ame**2 + p1p2)**2)
     &    - (C2apx35x12pm1p*1D0*
     &     (kp3*kp4*p1p4*p2p3 + kp3*kp4*p1p3*p2p4 + 
     &       ampi**2*kp4*p1p4*p2p4 + kp2*kp3*p1p4*p3p4 + 
     &       kp1*kp3*p2p4*p3p4 - kp3*p1p4*p2p4*p3p4*1D0 - 
     &       ame**2*ampi**4*kp4*2D0 + 
     &       ame**2*ampi**2*kp4**2*2D0 - 
     &       ampi**4*kp4*p1p2*2D0 + ampi**2*kp4**2*p1p2*2D0 - 
     &       ampi**2*kp2*kp4*p1p4*2D0 - 
     &       ampi**2*kp4*p1p4*p2p3*2D0 - 
     &       ampi**2*kp1*kp4*p2p4*2D0 - 
     &       ampi**2*kp4*p1p3*p2p4*2D0 + 
     &       ame**2*ampi**2*kp3*p3p4*2D0 + 
     &       ame**2*ampi**2*kp4*p3p4*2D0 - 
     &       ame**2*kp3*kp4*p3p4*2D0 + 
     &       ampi**2*kp3*p1p2*p3p4*2D0 + 
     &       ampi**2*kp4*p1p2*p3p4*2D0 - 
     &       kp3*kp4*p1p2*p3p4*2D0 + kp3*p1p4*p2p3*p3p4*2D0 + 
     &       kp3*p1p3*p2p4*p3p4*2D0 - 
     &       ame**2*kp3*p3p4**2*2D0 - kp3*p1p2*p3p4**2*2D0 + 
     &       ampi**2*kp1*kp2*kp4*3D0 + 
     &       ame**2*kp3**2*kp4*3D0 + kp3**2*kp4*p1p2*3D0 + 
     &       ampi**2*kp2*kp4*p1p3*3D0 - 
     &       kp2*kp3*kp4*p1p3*3D0 + 
     &       ampi**2*kp1*kp4*p2p3*3D0 - 
     &       kp1*kp3*kp4*p2p3*3D0 + 
     &       ampi**2*kp4*p1p3*p2p3*3D0 + 
     &       kp1*kp2*kp3*p3p4*3D0 - kp3*p1p3*p2p3*p3p4*3D0 - 
     &       kp3*kp4*p1p3*p2p3*6D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) - 
     &  (C2apx45x12pm1p*1D0*
     &     (ampi**2*kp3*p1p3*p2p3 + kp3*kp4*p1p4*p2p3 + 
     &       kp3*kp4*p1p3*p2p4 + kp2*kp4*p1p3*p3p4 + 
     &       kp1*kp4*p2p3*p3p4 - kp4*p1p3*p2p3*p3p4*1D0 - 
     &       ame**2*ampi**4*kp3*2D0 + 
     &       ame**2*ampi**2*kp3**2*2D0 - 
     &       ampi**4*kp3*p1p2*2D0 + ampi**2*kp3**2*p1p2*2D0 - 
     &       ampi**2*kp2*kp3*p1p3*2D0 - 
     &       ampi**2*kp1*kp3*p2p3*2D0 - 
     &       ampi**2*kp3*p1p4*p2p3*2D0 - 
     &       ampi**2*kp3*p1p3*p2p4*2D0 + 
     &       ame**2*ampi**2*kp3*p3p4*2D0 + 
     &       ame**2*ampi**2*kp4*p3p4*2D0 - 
     &       ame**2*kp3*kp4*p3p4*2D0 + 
     &       ampi**2*kp3*p1p2*p3p4*2D0 + 
     &       ampi**2*kp4*p1p2*p3p4*2D0 - 
     &       kp3*kp4*p1p2*p3p4*2D0 + kp4*p1p4*p2p3*p3p4*2D0 + 
     &       kp4*p1p3*p2p4*p3p4*2D0 - 
     &       ame**2*kp4*p3p4**2*2D0 - kp4*p1p2*p3p4**2*2D0 + 
     &       ampi**2*kp1*kp2*kp3*3D0 + 
     &       ame**2*kp3*kp4**2*3D0 + kp3*kp4**2*p1p2*3D0 + 
     &       ampi**2*kp2*kp3*p1p4*3D0 - 
     &       kp2*kp3*kp4*p1p4*3D0 + 
     &       ampi**2*kp1*kp3*p2p4*3D0 - 
     &       kp1*kp3*kp4*p2p4*3D0 + 
     &       ampi**2*kp3*p1p4*p2p4*3D0 + 
     &       kp1*kp2*kp4*p3p4*3D0 - kp4*p1p4*p2p4*p3p4*3D0 - 
     &       kp3*kp4*p1p4*p2p4*6D0))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D1ap0x12px35x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (ampi**2*kp4**2*(ame**2 + p1p2) + 
     &       kp3*p3p4*(kp2*p1p3 + kp1*(kp2 + p2p3) + 
     &          p2p4*(-(p1p4*1D0) + p1p3*2D0) + 
     &          (ame**2 + p1p2)*
     &           (-(kp3*1D0) + ampi**2*2D0 - p3p4*2D0) + 
     &          p2p3*(p1p4*2D0 - p1p3*3D0)) + 
     &       kp4*(kp3**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*kp2 - ampi**2*p2p4*1D0 + 
     &             p2p3*(-(kp3*1D0) + ampi**2*2D0)) + 
     &          ampi**2*(ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &             p1p2*(p3p4 - ampi**2*1D0)*2D0 + 
     &             p2p4*(p1p4 - p1p3*2D0) + 
     &             kp2*(-(p1p4*1D0) + p1p3*2D0) + 
     &             p2p3*(-(p1p4*2D0) + p1p3*3D0)) + 
     &          kp3*(p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p1p3*(p2p4 - p2p3*6D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D2ap0x12px35x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (ampi**2*kp4**2*x12 + 
     &       kp3*p3p4*(ame**2*(ampi**2 - p3p4*1D0)*2D0 + 
     &          p1p2*(ampi**2 - p3p4*1D0)*2D0 + 
     &          p1p4*(kp2 - p2p4*1D0 + p2p3*2D0) + 
     &          kp1*(p2p4 + kp2*3D0) + 
     &          p1p3*(p2p4*2D0 - p2p3*3D0)) + 
     &       kp4*(kp3**2*(ame**2 + p1p2)*3D0 + 
     &          kp1*(-(ampi**2*p2p4*2D0) + ampi**2*kp2*3D0 + 
     &             p2p3*(ampi**2 - kp3*1D0)*3D0) + 
     &          ampi**2*(ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &             p1p2*(p3p4 - ampi**2*1D0)*2D0 + 
     &             p2p4*(p1p4 - p1p3*2D0) + 
     &             (kp2 + p2p3)*(-(p1p4*2D0) + p1p3*3D0)) + 
     &          kp3*(p1p4*p2p3 - (ame**2 + p1p2)*p3p4*2D0 - 
     &             kp2*p1p3*3D0 + p1p3*(p2p4 - p2p3*6D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D1ap0x12px45x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*p3p4*(kp1*(kp2 + p2p4) + ampi**2*x12 - 
     &          p1p3*p2p3*1D0 + p1p3*p2p4*2D0 - 
     &          (ame**2 + p1p2)*1D0*(kp4 + p3p4*2D0) + 
     &          p1p4*(kp2 + p2p3*2D0 - p2p4*3D0)) + 
     &       kp3*(kp4**2*(ame**2 + p1p2) + 
     &          kp1*(ampi**2*(kp2 - p2p3*1D0) + 
     &             p2p4*(-(kp4*1D0) + ampi**2*2D0)) + 
     &          ampi**2*(ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &             p1p2*(p3p4 - ampi**2*1D0)*2D0 - 
     &             1D0*(kp2 - p2p3*1D0)*(p1p3 - p1p4*2D0) + 
     &             p2p4*(-(p1p3*2D0) + p1p4*3D0)) + 
     &          kp4*(p1p3*p2p4 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*6D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (D2ap0x12px45x34m1ppp*(m1**2 + p3p4*4D0)*
     &     (ampi**2*kp3**2*x12 + 
     &       kp4*p3p4*(kp2*p1p3 + x12*(ampi**2 - p3p4*1D0) - 
     &          p2p3*1D0*(p1p3 - p1p4*2D0) + 
     &          kp1*(p2p3 + kp2*3D0) + 
     &          p2p4*(p1p3*2D0 - p1p4*3D0)) + 
     &       kp3*(kp4**2*(ame**2 + p1p2)*3D0 + 
     &          kp1*(p2p4*(ampi**2 - kp4*1D0)*3D0 + 
     &             ampi**2*(-(p2p3*2D0) + kp2*3D0)) + 
     &          ampi**2*(ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &             p1p2*(p3p4 - ampi**2*1D0)*2D0 + 
     &             p2p3*(p1p3 - p1p4*2D0) + 
     &             kp2*(-(p1p3*2D0) + p1p4*3D0) + 
     &             p2p4*(-(p1p3*2D0) + p1p4*3D0)) + 
     &          kp4*(p1p3*p2p4 - (ame**2 + p1p2)*p3p4*2D0 + 
     &             p1p4*(p2p3 - kp2*3D0 - p2p4*6D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C12apx12x45m1pp*
     &     (-(kp3*1D0*(p1p4 + p2p4 - p3p4*2D0)*
     &          (ame**2*ampi**2 + ampi**2*p1p2 + kp2*p1p3 + 
     &            p1p3*p2p4 - 
     &            (ame**2 + p1p2)*(kp3 + p3p4)*1D0 + 
     &            p2p3*(kp1 + p1p4 - p1p3*2D0))) + 
     &       kp4*(-(1D0*(-(p1p3*1D0) - p2p3*1D0 + 
     &               ampi**2*2D0)*
     &             (p1p3*p2p4 - kp2*p1p3*1D0 + 
     &               ame**2*(ampi**2 - p3p4*1D0) + 
     &               p1p2*(ampi**2 - p3p4*1D0) + 
     &               p2p3*(p1p4 - kp1*1D0 - p1p3*2D0))) + 
     &          kp3*((ame**2 + p1p2)*(p1p3 + ampi**2*2D0) + 
     &             p2p3*(ame**2 + p1p2 - p1p3*8D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C12apx12x35m1pp*
     &     (-(kp4*1D0*(p1p3 + p2p3 - p3p4*2D0)*
     &          (ame**2*ampi**2 + ampi**2*p1p2 + 
     &            (kp1 + p1p3)*p2p4 - 
     &            (ame**2 + p1p2)*(kp4 + p3p4)*1D0 + 
     &            p1p4*(kp2 + p2p3 - p2p4*2D0))) + 
     &       kp3*(-(1D0*(-(p1p4*1D0) - p2p4*1D0 + 
     &               ampi**2*2D0)*
     &             (p2p4*(p1p3 - kp1*1D0) + 
     &               ame**2*(ampi**2 - p3p4*1D0) + 
     &               p1p2*(ampi**2 - p3p4*1D0) + 
     &               p1p4*(p2p3 - kp2*1D0 - p2p4*2D0))) + 
     &          kp4*((ame**2 + p1p2)*(p1p4 + ampi**2*2D0) + 
     &             p2p4*(ame**2 + p1p2 - p1p4*8D0)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C1ax45x12pm1pp*(-(kp3*1D0*
     &          (ame**2*ampi**2 + ampi**2*p1p2 + kp2*p1p3 + 
     &            p1p3*p2p4 - 
     &            (ame**2 + p1p2)*(kp3 + p3p4)*1D0 + 
     &            p2p3*(kp1 + p1p4 - p1p3*2D0))*
     &          (p1p4 + p2p4 - p3p4*3D0)) + 
     &       kp4*(-(1D0*(p1p3*p2p4 - kp2*p1p3*1D0 + 
     &               ame**2*(ampi**2 - p3p4*1D0) + 
     &               p1p2*(ampi**2 - p3p4*1D0) + 
     &               p2p3*(p1p4 - kp1*1D0 - p1p3*2D0))*
     &             (-(p1p3*1D0) - p2p3*1D0 + ampi**2*3D0)) + 
     &          kp3*((ame**2 + p1p2)*(p1p3 + ampi**2*3D0) + 
     &             p2p3*(ame**2 + p1p2 - p1p3*1.2D1)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2) + 
     &  (C1ax35x12pm1pp*(-(kp4*1D0*
     &          (ame**2*ampi**2 + ampi**2*p1p2 + 
     &            (kp1 + p1p3)*p2p4 - 
     &            (ame**2 + p1p2)*(kp4 + p3p4)*1D0 + 
     &            p1p4*(kp2 + p2p3 - p2p4*2D0))*
     &          (p1p3 + p2p3 - p3p4*3D0)) + 
     &       kp3*(-(1D0*(p2p4*(p1p3 - kp1*1D0) + 
     &               ame**2*(ampi**2 - p3p4*1D0) + 
     &               p1p2*(ampi**2 - p3p4*1D0) + 
     &               p1p4*(p2p3 - kp2*1D0 - p2p4*2D0))*
     &             (-(p1p4*1D0) - p2p4*1D0 + ampi**2*3D0)) + 
     &          kp4*((ame**2 + p1p2)*(p1p4 + ampi**2*3D0) + 
     &             p2p4*(ame**2 + p1p2 - p1p4*1.2D1)))))/
     &   (kp3*kp4*(ame**2 + p1p2)**2)
