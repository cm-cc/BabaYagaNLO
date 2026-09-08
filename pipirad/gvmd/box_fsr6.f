        elmat2 = (D1ax15x34px24epem1m2p*(ampi**2 - p2p4*1D0)*2D0*
     &     (kp3**2*(ampi**2*(ame**2 + p1p2) - 
     &          p1p4*1D0*(ame**2 + p1p2 + p2p4 - p2p3*1D0) - 
     &          ame**2*kp4*2D0) + 
     &       kp4*(kp2*(p1p3*(p1p3 + p1p4 + p3p4) - 
     &             ampi**2*p1p4*1D0) + 
     &          kp4*(ame**2*ampi**2 + 
     &             p1p2*(ampi**2 - p1p3*1D0) - 
     &             p1p3*1D0*(ame**2 + p2p3 - p2p4*1D0)) + 
     &          p1p3*(p1p3*p2p3 - p1p4*p2p4*1D0)*2D0) + 
     &       kp3*(-(ampi**2*kp2*p1p3*1D0) - 
     &          ame**2*kp4**2*2D0 + 
     &          p1p4*(kp2*(p1p3 + p1p4 + p3p4) - 
     &             p1p3*p2p3*2D0 + p1p4*p2p4*2D0) + 
     &          kp4*(p1p3*
     &              (ame**2 + p2p4 - p1p2*1D0 - p2p3*1D0) + 
     &             p1p4*(ame**2 + p2p3 - p1p2*1D0 - 
     &                p2p4*1D0) - ame**2*(p2p3 + p2p4)*2D0 - 
     &             (ame**2 + p1p2)*p3p4*2D0)) + 
     &       kp1*(kp2*(kp3 + kp4)*(ampi**2 + p3p4) - 
     &          kp4**2*p2p3*1D0 + 
     &          kp3*(p1p4*(p2p3 + p2p4) - p2p3*p3p4*1D0 + 
     &             p2p4*(ampi**2 - kp3*1D0)) + 
     &          kp4*(p2p3*(ampi**2 + p1p3 - kp3*3D0) + 
     &             p2p4*(p1p3 - p3p4*1D0 - kp3*3D0)))))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) - 
     &  (D0ax15x34px24epem1m2p*p2p4*2D0*
     &     (kp3**2*(p1p4*(ame**2 + p1p2 + p2p4 - p2p3*1D0) + 
     &          ame**2*kp4*2D0) - 
     &       kp4*p1p3*1D0*
     &        (kp2*(ampi**2 + p1p3 - p1p4*1D0 - p3p4*1D0) + 
     &          (p2p3 - p2p4*1D0)*
     &           (kp4 + p1p3*2D0 - p1p4*2D0) + 
     &          ame**2*(kp4 - ampi**2*2D0 + p3p4*2D0) + 
     &          p1p2*(kp4 - ampi**2*2D0 + p3p4*2D0)) + 
     &       kp3*(-(ame**2*kp4**2*2D0) + 
     &          kp4*(p1p3*
     &              (p1p2 + p2p3 - ame**2*1D0 - p2p4*1D0) + 
     &             p1p4*(ame**2 + p2p3 - p1p2*1D0 - 
     &                p2p4*1D0) + ame**2*(p2p3 - p2p4*1D0)*2D0
     &             ) + p1p4*
     &           (p3p4*x12 + 
     &             kp2*(ampi**2 + p1p4 - p1p3*1D0 - 
     &                p3p4*1D0) - ame**2*ampi**2*2D0 - 
     &             ampi**2*p1p2*2D0 + 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0))
     &        + kp1*(kp3**2*p2p4 + 
     &          kp2*(kp3 - kp4*1D0)*(ampi**2 - p3p4*1D0) - 
     &          kp4*1D0*(kp4*p2p3 + 
     &             p2p4*(p3p4 - ampi**2*1D0) + 
     &             p1p3*(p2p3 - p2p4*1D0)) + 
     &          kp3*(p2p4*(p1p4 - kp4*3D0) - 
     &             p2p3*1D0*
     &              (ampi**2 + p1p4 - p3p4*1D0 - kp4*3D0)))))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D3ax15x34px24epem1m2p*
     &     (ampi**2*kp3**2*
     &        (ampi**2*(ame**2 + p1p2) - 
     &          (ame**2 + p1p2)*p2p4*1D0 - 
     &          p1p4*1D0*(ame**2 + p1p2 + p2p4 - p2p3*1D0) - 
     &          ame**2*kp4*2D0) + 
     &       kp4*(kp2*(ampi**2*p1p3**2 + 
     &             ampi**2*p1p4*(p2p4 - ampi**2*1D0) + 
     &             p1p3*(ampi**2*p2p4 + 
     &                p1p4*(ampi**2 - p2p4*2D0) + 
     &                p3p4*(ampi**2 - p2p4*2D0))) + 
     &          kp4*(ame**2*ampi**2*(ampi**2 - p2p4*1D0) - 
     &             p1p3*1D0*(ame**2 + p2p3 - p2p4*1D0)*
     &              (ampi**2 - p2p4*2D0) + 
     &             p1p2*(ampi**4 - ampi**2*p2p4*1D0 - 
     &                p1p3*1D0*(ampi**2 - p2p4*2D0))) - 
     &          p1p3*2D0*
     &           (p1p3*(p2p4**2 - ampi**2*p2p3*1D0) + 
     &             p2p4*(ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0) + 
     &                p1p4*(ampi**2 + p2p3 - p2p4*2D0)))) + 
     &       kp3*(-(ame**2*kp4**2*2D0*(ampi**2 - p2p4*2D0)) + 
     &          kp2*(ampi**2*p1p3*
     &              (p1p4 + p2p4 - ampi**2*1D0) + 
     &             p1p4*(ampi**2*(p3p4 - p2p4*1D0) + 
     &                p1p4*(ampi**2 - p2p4*2D0))) + 
     &          kp4*(ampi**2*p1p3*
     &              (ame**2 + p2p4 - p1p2*1D0 - p2p3*1D0) - 
     &             (ame**2 + p1p2)*p3p4*(ampi**2 - p2p4*1D0)*
     &              2D0 + 
     &             p1p4*(ame**2 + p2p3 - p1p2*1D0 - p2p4*1D0)*
     &              (ampi**2 - p2p4*2D0) - 
     &             ame**2*2D0*
     &              (ampi**2*p2p3 + p2p4*(ampi**2 - p2p4*2D0))
     &             ) + p1p4*2D0*
     &           (p1p3*(p2p4**2 - ampi**2*p2p3*1D0) + 
     &             p2p4*(ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0) + 
     &                p1p4*(ampi**2 + p2p3 - p2p4*2D0)))) + 
     &       kp1*(-(kp4**2*p2p3*1D0*(ampi**2 - p2p4*2D0)) + 
     &          kp3*(ampi**2*
     &              (-(p2p3*p3p4*1D0) + 
     &                p2p4*
     &                 (ampi**2 + p2p3 - kp3*1D0 - p2p4*1D0))
     &              + p1p4*
     &              (ampi**2*p2p3 + p2p4*(ampi**2 - p2p4*2D0))
     &             ) + kp2*
     &           (ampi**2*kp3*(ampi**2 + p3p4 - p2p4*2D0) + 
     &             kp4*(ampi**4 + p3p4*(ampi**2 - p2p4*2D0)))
     &           + kp4*(ampi**2*p2p3*
     &              (ampi**2 + p1p3 - p2p4*1D0 - kp3*3D0) - 
     &             p2p4*1D0*
     &              (ampi**2*p2p4 + 
     &                p3p4*(ampi**2 - p2p4*2D0) - 
     &                p1p3*1D0*(ampi**2 - p2p4*2D0) + 
     &                kp3*(ampi**2 - p2p4*2D0)*3D0)))))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) + 
     &  (C0appx34m1pm2*2.5D-1*
     &     (kp3**2*(ampi**2*(ame**2 + p1p2) - kp1*p2p4*2D0 - 
     &          p1p4*(ame**2 + p1p2 + p2p4 - p2p3*1D0)*2D0 - 
     &          ame**2*kp4*4D0) + 
     &       kp4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          ampi**2*kp2*p1p3 + p1p3*p3p4*x12 - 
     &          ampi**2*kp2*p1p4*1D0 - 
     &          ame**2*ampi**2*p1p3*2D0 - 
     &          ampi**2*p1p2*p1p3*2D0 + kp2*p1p3**2*2D0 - 
     &          p1p3*p1p4*p2p3*2D0 - p1p3**2*p2p4*2D0 + 
     &          kp1*(-(ampi**2*p2p4*1D0) + ampi**2*kp2*2D0 + 
     &             p2p3*(ampi**2 + p1p3*2D0)) + 
     &          p1p3**2*p2p3*4D0) + 
     &       kp3*(ampi**2*kp1*p2p3 + ampi**2*kp1*p2p4 + 
     &          ame**2*ampi**2*p1p4*2D0 + 
     &          ampi**2*p1p2*p1p4*2D0 + kp1*p1p4*p2p3*2D0 + 
     &          p1p4**2*p2p3*2D0 + p1p3*p1p4*p2p4*2D0 - 
     &          ((ame**2 + p1p2)*p1p4 + kp1*p2p3)*p3p4*2D0 + 
     &          kp2*(-(ampi**2*p1p4*1D0) + 
     &             (kp1 + p1p4)*p3p4*2D0 - 
     &             p1p3*1D0*(ampi**2 - p1p4*2D0)) + 
     &          kp4*2D0*(-((ame**2 + p1p2)*p3p4*1D0) + 
     &             p1p3*(ame**2 + p2p4 - p1p2*1D0 - 
     &                p2p3*1D0) - 
     &             p2p3*1D0*(ame**2*2D0 + kp1*3D0)) - 
     &          p1p3*p1p4*p2p3*4D0)))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) - 
     &  (C1apx34ppm1m2*2.5D-1*
     &     (kp3**2*(ampi**2*(ame**2 + p1p2) - kp1*p2p4*2D0 - 
     &          p1p4*(ame**2 + p1p2 + p2p4 - p2p3*1D0)*2D0 - 
     &          ame**2*kp4*4D0) + 
     &       kp4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &          ampi**2*kp2*p1p3 + p1p3*p3p4*x12 - 
     &          ampi**2*kp2*p1p4*1D0 - 
     &          ame**2*ampi**2*p1p3*2D0 - 
     &          ampi**2*p1p2*p1p3*2D0 + kp2*p1p3**2*2D0 - 
     &          p1p3*p1p4*p2p3*2D0 - p1p3**2*p2p4*2D0 + 
     &          kp1*(-(ampi**2*p2p4*1D0) + ampi**2*kp2*2D0 + 
     &             p2p3*(ampi**2 + p1p3*2D0)) + 
     &          p1p3**2*p2p3*4D0) + 
     &       kp3*(ampi**2*kp1*p2p3 + ampi**2*kp1*p2p4 + 
     &          ame**2*ampi**2*p1p4*2D0 + 
     &          ampi**2*p1p2*p1p4*2D0 + kp1*p1p4*p2p3*2D0 + 
     &          p1p4**2*p2p3*2D0 + p1p3*p1p4*p2p4*2D0 - 
     &          ((ame**2 + p1p2)*p1p4 + kp1*p2p3)*p3p4*2D0 + 
     &          kp2*(-(ampi**2*p1p4*1D0) + 
     &             (kp1 + p1p4)*p3p4*2D0 - 
     &             p1p3*1D0*(ampi**2 - p1p4*2D0)) + 
     &          kp4*2D0*(-((ame**2 + p1p2)*p3p4*1D0) + 
     &             p1p3*(ame**2 + p2p4 - p1p2*1D0 - 
     &                p2p3*1D0) - 
     &             p2p3*1D0*(ame**2*2D0 + kp1*3D0)) - 
     &          p1p3*p1p4*p2p3*4D0)))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) + 
     &  (C1apx34ppm2m1*2.5D-1*
     &     (ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*(-(ampi**2*kp2*p1p4*1D0) + 
     &          p1p3**2*p2p4*2D0 + 
     &          kp4*(ame**2*ampi**2 - 
     &             p1p3*(ame**2 + p2p3 - p2p4*1D0)*2D0 + 
     &             p1p2*(ampi**2 - p1p3*2D0)) + 
     &          p1p3*(kp2*
     &              (-(ampi**2*1D0) + p1p4*2D0 + p3p4*2D0) + 
     &             2D0*(ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0) + 
     &                p1p4*(p2p3 - p2p4*2D0)))) + 
     &       kp3*(-(ampi**2*kp2*p1p3*1D0) + 
     &          kp4*2D0*(-((ame**2 + p1p2)*p3p4*1D0) + 
     &             p1p4*(ame**2 + p2p3 - p1p2*1D0 - 
     &                p2p4*1D0) - ame**2*p2p4*2D0) + 
     &          p1p4*(p3p4*x12 + kp2*(ampi**2 + p1p4*2D0) - 
     &             2D0*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &                p1p3*p2p4 + p1p4*(p2p3 - p2p4*2D0))) - 
     &          ame**2*kp4**2*4D0) + 
     &       kp1*(kp2*(ampi**2*kp3 + kp4*p3p4)*2D0 + 
     &          kp4*(p2p3*(ampi**2 - kp4*2D0) + 
     &             p2p4*(ampi**2 + p1p3*2D0 - p3p4*2D0)) + 
     &          kp3*(-(ampi**2*p2p3*1D0) + 
     &             p2p4*(ampi**2 + p1p4*2D0 - kp4*6D0)))))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) - 
     &  (C0aex34x15em2m1*5.D-1*
     &     (ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*(-(ampi**2*kp2*p1p4*1D0) + 
     &          p1p3**2*p2p4*2D0 + 
     &          kp4*(ame**2*ampi**2 - 
     &             p1p3*(ame**2 + p2p3 - p2p4*1D0)*2D0 + 
     &             p1p2*(ampi**2 - p1p3*2D0)) + 
     &          p1p3*(kp2*
     &              (-(ampi**2*1D0) + p1p4*2D0 + p3p4*2D0) + 
     &             2D0*(ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0) + 
     &                p1p4*(p2p3 - p2p4*2D0)))) + 
     &       kp3*(-(ampi**2*kp2*p1p3*1D0) + 
     &          kp4*2D0*(-((ame**2 + p1p2)*p3p4*1D0) + 
     &             p1p4*(ame**2 + p2p3 - p1p2*1D0 - 
     &                p2p4*1D0) - ame**2*p2p4*2D0) + 
     &          p1p4*(p3p4*x12 + kp2*(ampi**2 + p1p4*2D0) - 
     &             2D0*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &                p1p3*p2p4 + p1p4*(p2p3 - p2p4*2D0))) - 
     &          ame**2*kp4**2*4D0) + 
     &       kp1*(kp2*(ampi**2*kp3 + kp4*p3p4)*2D0 + 
     &          kp4*(p2p3*(ampi**2 - kp4*2D0) + 
     &             p2p4*(ampi**2 + p1p3*2D0 - p3p4*2D0)) + 
     &          kp3*(-(ampi**2*p2p3*1D0) + 
     &             p2p4*(ampi**2 + p1p4*2D0 - kp4*6D0)))))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D1aex24px34px15m2epm1*
     &     (ame**2*ampi**2*kp3**2*p1p4 + 
     &       ame**2*kp3**2*p1p4**2 + 
     &       ame**2*ampi**2*kp3**2*p2p3 + 
     &       ame**2*ampi**2*kp4**2*p2p3 + 
     &       ame**2*kp3*kp4*p1p4*p2p3 + 
     &       ame**2*ampi**2*kp3**2*p2p4 + 
     &       ampi**2*kp3**2*p1p2*p2p4 + 
     &       ampi**2*kp4**2*p1p2*p2p4 + 
     &       ame**2*kp3*kp4*p1p4*p2p4 + 
     &       ame**2*kp3*kp4*p2p4**2 - 
     &       ame**2*ampi**2*kp3*kp4*p1p3*1D0 - 
     &       ame**2*kp4**2*p1p3**2*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*p2p3*1D0 - 
     &       ame**2*kp4**2*p1p3*p2p3*1D0 - 
     &       ame**2*kp4**2*p1p3*p2p4*1D0 - 
     &       ame**2*kp4**2*p2p3*p2p4*1D0 - 
     &       ame**2*kp3**2*p2p4**2*1D0 - 
     &       ame**4*ampi**2*kp3*kp4*2D0 + 
     &       ame**2*ampi**2*kp3*kp4*p1p2*2D0 + 
     &       ame**2*ampi**2*kp4*p1p3**2*2D0 - 
     &       ame**2*ampi**2*kp3*p1p3*p1p4*2D0 - 
     &       ame**2*kp3*kp4**2*p2p3*2D0 - 
     &       ame**2*kp3**2*kp4*p2p4*2D0 + 
     &       ame**2*kp3*kp4*p1p3*p2p4*2D0 - 
     &       ampi**2*kp4*p1p2*p1p3*p2p4*2D0 - 
     &       kp3*kp4*p1p2*p1p3*p2p4*2D0 - 
     &       ame**2*kp4*p1p3**2*p2p4*2D0 - 
     &       ame**2*kp3**2*p1p4*p2p4*2D0 + 
     &       ampi**2*kp3*p1p2*p1p4*p2p4*2D0 - 
     &       kp3**2*p1p2*p1p4*p2p4*2D0 + 
     &       ame**2*kp3*p1p3*p1p4*p2p4*2D0 + 
     &       ame**2*kp4*p1p3*p1p4*p2p4*2D0 - 
     &       ame**2*kp3*p1p4**2*p2p4*2D0 + 
     &       ame**2*kp4*p1p3*p2p3*p2p4*2D0 - 
     &       kp3*kp4*p1p3*p2p3*p2p4*2D0 - 
     &       ame**2*kp3*p1p4*p2p3*p2p4*2D0 + 
     &       kp3**2*p1p4*p2p3*p2p4*2D0 - 
     &       kp4*p1p3*p1p4*p2p3*p2p4*2D0 + 
     &       kp3*p1p4**2*p2p3*p2p4*2D0 - 
     &       ame**2*kp4*p1p3*p2p4**2*2D0 + 
     &       kp3*kp4*p1p3*p2p4**2*2D0 - 
     &       kp4*p1p3**2*p2p4**2*2D0 + 
     &       ame**2*kp3*p1p4*p2p4**2*2D0 - 
     &       kp3**2*p1p4*p2p4**2*2D0 + 
     &       kp3*p1p3*p1p4*p2p4**2*2D0 + 
     &       p3p4*(-(ame**2*kp3**2*(p1p4 + p2p4)*1D0) + 
     &          kp4*p1p3*(p1p2*p2p4 - ame**2*p1p3*1D0)*2D0 + 
     &          kp3*(-(p1p2*(kp4 + p1p4)*p2p4*2D0) + 
     &             ame**2*
     &              (p1p3*p1p4*2D0 + 
     &                kp4*
     &                 (p1p3 - p2p3*1D0 + ame**2*2D0 - 
     &                   p1p2*2D0)))) + 
     &       kp2*(kp4*(ame**2*ampi**4 - 
     &             ampi**2*(ame**2 + p1p4)*p2p4*1D0 + 
     &             ame**2*p3p4*(p2p4 - ampi**2*1D0) + 
     &             p1p3**2*p2p4*2D0 + 
     &             ame**2*kp3*2D0*
     &              (p3p4 - ampi**2*1D0 + p2p4*2D0) + 
     &             p1p3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &                p2p4*(ampi**2 + ame**2*2D0))) + 
     &          kp3*(ame**2*p3p4*(p3p4 - ampi**2*1D0) + 
     &             p2p4*(ame**2*p3p4 - 
     &                ampi**2*(ame**2 + p1p3)*1D0) + 
     &             p1p4*(ame**2*(p3p4 - ampi**2*1D0) + 
     &                p2p4*
     &                 (-(ampi**2*1D0) + ame**2*2D0 + 
     &                   p1p3*2D0 + p3p4*2D0)))) - 
     &       ame**2*kp3*kp4*p2p3*p2p4*3D0 + 
     &       ame**4*kp3*kp4*p2p4*4D0 - 
     &       ame**2*kp3*kp4*p1p2*p2p4*4D0 + 
     &       kp4*p1p3**2*p2p3*p2p4*4D0 - 
     &       kp3*p1p3*p1p4*p2p3*p2p4*4D0 + 
     &       kp1*(kp3**2*(ame**2*ampi**2 - p2p4**2*2D0) + 
     &          kp4*(ame**2*ampi**4 - 
     &             ame**2*ampi**2*kp4*1D0 + 
     &             ame**2*p3p4*(p2p4 - ampi**2*1D0) - 
     &             ampi**2*p2p4*1D0*(p2p4 + x25 - p2p3*1D0) + 
     &             p1p3*(ame**2*(ampi**2 - p3p4*1D0) + 
     &                p2p4*(p2p3 - ame**2*1D0)*2D0)) + 
     &          kp3*(ame**2*p3p4**2 + 
     &             ampi**2*p2p4*(p2p3 + p2p4 - ame**2*1D0) + 
     &             p1p4*(ame**2*(ampi**2 - p3p4*1D0) + 
     &                p2p4*(p2p3 - ame**2*1D0)*2D0) + 
     &             p3p4*(-(ame**2*ampi**2*1D0) + 
     &                p2p4*(ame**2 + kp2*2D0 - p2p3*2D0)) + 
     &             kp4*(ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &                p2p4*(ame**2*4D0 - p2p3*6D0))))))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2))
