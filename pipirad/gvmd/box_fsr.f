        elmat2 = (C1apx45x12m1pm2*ep3*5.D-1*
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
     &   (kp3*kp4**2*(ame**2 + p1p2)) - 
     &  (C1apx45x12m2pm1*ep3*5.D-1*
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
     &   (kp3*kp4**2*(ame**2 + p1p2)) - 
     &  (C1apx35x12m1pm2*ep3*5.D-1*
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
     &   (kp3**2*kp4*(ame**2 + p1p2)) + 
     &  (C1apx35x12m2pm1*ep3*5.D-1*
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
     &   (kp3**2*kp4*(ame**2 + p1p2)) + 
     &  (C0ax15ex34m1em2*ep2*2D0*
     &     (kp4*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)) + 
     &  (D0aex34px13x25pem1m2p*ep2*p1p3*2D0*
     &     (kp1*(kp4*p2p3*
     &           (p2p4 + p3p4 - ampi**2*1D0 - p2p3*1D0) + 
     &          kp2*(kp3 - kp4*1D0)*(ampi**2 - p3p4*1D0) + 
     &          kp3*p2p4*
     &           (ampi**2 + p2p4 - p2p3*1D0 - p3p4*1D0)) + 
     &       kp3**2*(kp2*p1p4 + 
     &          p2p4*(ame**2 + p1p2 + p1p4 - p1p3*1D0) + 
     &          ame**2*kp4*2D0) - 
     &       kp4*1D0*(kp2*
     &           (kp4*p1p3 + p1p3*p2p3 - 
     &             p1p4*1D0*(ampi**2 + p2p3 - p3p4*1D0)) + 
     &          p2p3*(p3p4*x12 + 
     &             kp4*(ame**2 + p1p2 + p1p3 - p1p4*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0))
     &        + kp3*(kp2*
     &           (p1p4*p2p4 - 
     &             p1p3*1D0*(ampi**2 + p2p4 - p3p4*1D0)) - 
     &          ame**2*kp4**2*2D0 + 
     &          p2p4*(p3p4*x12 - ame**2*ampi**2*2D0 - 
     &             ampi**2*p1p2*2D0 + 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0) + 
     &          kp4*(-(1D0*(ame**2 - p1p2*1D0)*
     &                (p2p3 - p2p4*1D0)) + 
     &             p1p3*
     &              (p2p3 + p2p4 + ame**2*2D0 + kp2*3D0) - 
     &             p1p4*1D0*
     &              (p2p3 + p2p4 + ame**2*2D0 + kp2*3D0)))))/
     &   (kp2*kp3*kp4*(ame**2 + p1p2)) - 
     &  (D0aex34px14x25pem1m2p*ep2*p1p4*2D0*
     &     (kp1*(kp4*p2p3*
     &           (p2p4 + p3p4 - ampi**2*1D0 - p2p3*1D0) + 
     &          kp2*(kp3 - kp4*1D0)*(ampi**2 - p3p4*1D0) + 
     &          kp3*p2p4*
     &           (ampi**2 + p2p4 - p2p3*1D0 - p3p4*1D0)) + 
     &       kp3**2*(kp2*p1p4 + 
     &          p2p4*(ame**2 + p1p2 + p1p4 - p1p3*1D0) + 
     &          ame**2*kp4*2D0) - 
     &       kp4*1D0*(kp2*
     &           (kp4*p1p3 + p1p3*p2p3 - 
     &             p1p4*1D0*(ampi**2 + p2p3 - p3p4*1D0)) + 
     &          p2p3*(p3p4*x12 + 
     &             kp4*(ame**2 + p1p2 + p1p3 - p1p4*1D0) - 
     &             ampi**2*(ame**2 + p1p2)*2D0 + 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0))
     &        + kp3*(kp2*
     &           (p1p4*p2p4 - 
     &             p1p3*1D0*(ampi**2 + p2p4 - p3p4*1D0)) - 
     &          ame**2*kp4**2*2D0 + 
     &          p2p4*(p3p4*x12 - ame**2*ampi**2*2D0 - 
     &             ampi**2*p1p2*2D0 + 
     &             (p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0)*2D0) + 
     &          kp4*(-(1D0*(ame**2 - p1p2*1D0)*
     &                (p2p3 - p2p4*1D0)) + 
     &             p1p3*
     &              (p2p3 + p2p4 + ame**2*2D0 + kp2*3D0) - 
     &             p1p4*1D0*
     &              (p2p3 + p2p4 + ame**2*2D0 + kp2*3D0)))))/
     &   (kp2*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D2aex34px13x25pem1m2p*ep2*(ampi**2 - p1p3*1D0)*2D0*
     &     (kp1*(kp2*(kp3 + kp4)*(ampi**2 + p3p4) + 
     &          kp4*(p2p3*(p2p3 + p2p4 + p3p4) - 
     &             ampi**2*p2p4*1D0) + 
     &          kp3*(p2p4*(p2p4 + p3p4) + 
     &             p2p3*(p2p4 - ampi**2*1D0))) + 
     &       kp3**2*(ampi**2*(ame**2 + p1p2) - kp2*p1p4*1D0 - 
     &          p2p4*1D0*(ame**2 + p1p2 + p1p4 - p1p3*1D0) - 
     &          ame**2*kp4*2D0) + 
     &       kp4*(kp4*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp2*p1p3*1D0 - 
     &             p2p3*1D0*(ame**2 + p1p2 + p1p3 - p1p4*1D0))
     &            + kp2*(p1p3*(ampi**2 + p2p3) + 
     &             p1p4*(p2p3 - p3p4*1D0)) + 
     &          p2p3*(p1p3*p2p3 - p1p4*p2p4*1D0)*2D0) - 
     &       kp3*1D0*(-(kp2*1D0*
     &             (p1p4*(ampi**2 + p2p4) + 
     &               p1p3*(p2p4 - p3p4*1D0))) + 
     &          ame**2*kp4**2*2D0 + 
     &          p2p4*(p1p3*p2p3 - p1p4*p2p4*1D0)*2D0 + 
     &          kp4*(p3p4*x12 - 
     &             (p2p3 + p2p4)*1D0*(ame**2 - p1p2*1D0) + 
     &             p1p4*(p2p4 - p2p3*1D0 + ame**2*2D0 + 
     &                kp2*3D0) + 
     &             p1p3*(p2p3 - p2p4*1D0 + ame**2*2D0 + 
     &                kp2*3D0)))))/
     &   (kp2*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D2aex34px14x25pem1m2p*ep2*(ampi**2 - p1p4*1D0)*2D0*
     &     (kp1*(kp2*(kp3 + kp4)*(ampi**2 + p3p4) + 
     &          kp4*(p2p3*(p2p3 + p2p4 + p3p4) - 
     &             ampi**2*p2p4*1D0) + 
     &          kp3*(p2p4*(p2p4 + p3p4) + 
     &             p2p3*(p2p4 - ampi**2*1D0))) + 
     &       kp3**2*(ampi**2*(ame**2 + p1p2) - kp2*p1p4*1D0 - 
     &          p2p4*1D0*(ame**2 + p1p2 + p1p4 - p1p3*1D0) - 
     &          ame**2*kp4*2D0) + 
     &       kp4*(kp4*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp2*p1p3*1D0 - 
     &             p2p3*1D0*(ame**2 + p1p2 + p1p3 - p1p4*1D0))
     &            + kp2*(p1p3*(ampi**2 + p2p3) + 
     &             p1p4*(p2p3 - p3p4*1D0)) + 
     &          p2p3*(p1p3*p2p3 - p1p4*p2p4*1D0)*2D0) - 
     &       kp3*1D0*(-(kp2*1D0*
     &             (p1p4*(ampi**2 + p2p4) + 
     &               p1p3*(p2p4 - p3p4*1D0))) + 
     &          ame**2*kp4**2*2D0 + 
     &          p2p4*(p1p3*p2p3 - p1p4*p2p4*1D0)*2D0 + 
     &          kp4*(p3p4*x12 - 
     &             (p2p3 + p2p4)*1D0*(ame**2 - p1p2*1D0) + 
     &             p1p4*(p2p4 - p2p3*1D0 + ame**2*2D0 + 
     &                kp2*3D0) + 
     &             p1p3*(p2p3 - p2p4*1D0 + ame**2*2D0 + 
     &                kp2*3D0)))))/
     &   (kp2*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D1ax15x34px23epem1m2p*ep2*(ampi**2 - p2p3*1D0)*2D0*
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
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D1ax15x34px24epem1m2p*ep2*(ampi**2 - p2p4*1D0)*2D0*
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
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D0ax15x34px23epem1m2p*ep2*p2p3*2D0*
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
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) - 
     &  (D0ax15x34px24epem1m2p*ep2*p2p4*2D0*
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
     &  (D3ax15x34px24epem1m2p*ep2*
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
     &  (D3aex34px14x25pem1m2p*ep2*
     &     (ampi**2*kp3**2*
     &        (ampi**2*(ame**2 + p1p2) - 
     &          (ame**2 + kp2 + p1p2)*p1p4*1D0 - 
     &          p2p4*1D0*(ame**2 + p1p2 + p1p4 - p1p3*1D0) - 
     &          ame**2*kp4*2D0) + 
     &       kp4*(kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &             ampi**2*kp2*p1p3*1D0 - 
     &             ame**2*ampi**2*p1p4*1D0 - 
     &             ampi**2*p1p2*p1p4*1D0 + 
     &             kp2*p1p3*p1p4*2D0 - 
     &             p2p3*1D0*(ame**2 + p1p2 + p1p3 - p1p4*1D0)*
     &              (ampi**2 - p1p4*2D0)) + 
     &          p2p3*2D0*
     &           (ampi**2*p1p3*p2p3 + 
     &             p1p4*((ame**2 + p1p2)*p3p4 - 
     &                ampi**2*(ame**2 + p1p2)*1D0 - 
     &                (ampi**2 + p1p3)*p2p4*1D0) - 
     &             p1p4**2*1D0*(p2p3 - p2p4*2D0)) + 
     &          kp2*(ampi**2*p1p3*
     &              (ampi**2 + p2p3 - p1p4*1D0) - 
     &             p1p4*1D0*
     &              (ampi**2*(p3p4 - p2p3*1D0) + 
     &                p1p4*(ampi**2 + p2p3*2D0 - p3p4*2D0))))
     &        + kp1*(kp2*
     &           (ampi**2*kp3*(ampi**2 + p3p4 - p1p4*2D0) + 
     &             kp4*(ampi**4 + p3p4*(ampi**2 - p1p4*2D0)))
     &           + kp3*(ampi**2*p2p3*
     &              (p1p4 + p2p4 - ampi**2*1D0) + 
     &             p2p4*(ampi**2*(p2p4 + p3p4) - 
     &                p1p4*1D0*(ampi**2 + p2p4*2D0))) + 
     &          kp4*(ampi**2*p2p3*(p2p3 + p2p4 + p3p4) - 
     &             ampi**4*p2p4*1D0 + 
     &             p1p4*(ampi**2*p2p4 + 
     &                p2p3*(ampi**2 - p2p4*2D0 - p3p4*2D0))))
     &        + kp3*(-(ame**2*kp4**2*2D0*
     &             (ampi**2 - p1p4*2D0)) + 
     &          p2p4*2D0*
     &           (-(ampi**2*p1p3*p2p3*1D0) + 
     &             p1p4*((ampi**2 + p1p3)*p2p4 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0)) + 
     &             p1p4**2*(p2p3 - p2p4*2D0)) + 
     &          kp2*(ampi**2*p1p4*(ampi**2 + p1p3 + p2p4) + 
     &             ampi**2*p1p3*(p2p4 - p3p4*1D0) - 
     &             p1p4**2*1D0*(ampi**2 + p2p4*2D0)) - 
     &          kp4*1D0*(ampi**2*
     &              (-((p2p3 + p2p4)*1D0*
     &                   (ame**2 - p1p2*1D0)) + 
     &                (ame**2 + p1p2)*p3p4*2D0) - 
     &             p1p4**2*2D0*
     &              (p2p4 + ame**2*2D0 + kp2*3D0) + 
     &             p1p4*(-(p3p4*x12*1D0) + 
     &                ame**2*(ampi**2 + p2p4)*2D0 + 
     &                p2p4*(ampi**2 - p1p2*2D0) + 
     &                ampi**2*(-(p2p3*1D0) + kp2*3D0)) + 
     &             p1p3*(-(p2p4*1D0*(ampi**2 - p1p4*2D0)) + 
     &                ampi**2*(p2p3 + ame**2*2D0 + kp2*3D0))))
     &       ))/(kp2*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D3aex34px13x25pem1m2p*ep2*
     &     (kp3**2*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &          ame**2*ampi**2*p1p3*1D0 - 
     &          ampi**2*p1p2*p1p3*1D0 - 
     &          ampi**2*kp2*p1p4*1D0 + kp2*p1p3*p1p4*2D0 - 
     &          p2p4*1D0*(ame**2 + p1p2 + p1p4 - p1p3*1D0)*
     &           (ampi**2 - p1p3*2D0) - 
     &          ame**2*kp4*2D0*(ampi**2 - p1p3*2D0)) + 
     &       kp4*(ampi**2*kp4*
     &           (ame**2*ampi**2 - 
     &             p1p3*(ame**2 + kp2 + p2p3)*1D0 + 
     &             p2p3*(p1p4 - ame**2*1D0) + 
     &             p1p2*(ampi**2 - p1p3*1D0 - p2p3*1D0)) + 
     &          p2p3*2D0*
     &           (-(ampi**2*p1p4*p2p4*1D0) + 
     &             p1p3*((ampi**2 + p1p4)*p2p3 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0)) + 
     &             p1p3**2*(p2p4 - p2p3*2D0)) + 
     &          kp2*(ampi**2*p1p3*(ampi**2 + p1p4 + p2p3) + 
     &             ampi**2*p1p4*(p2p3 - p3p4*1D0) - 
     &             p1p3**2*1D0*(ampi**2 + p2p3*2D0))) + 
     &       kp1*(kp2*(ampi**2*kp4*
     &              (ampi**2 + p3p4 - p1p3*2D0) + 
     &             kp3*(ampi**4 + p3p4*(ampi**2 - p1p3*2D0)))
     &           + kp3*(p2p3*
     &              (ampi**2*p1p3 - ampi**4*1D0 + 
     &                p2p4*(ampi**2 - p1p3*2D0)) + 
     &             p2p4*(ampi**2*(p1p3 + p2p4) + 
     &                p3p4*(ampi**2 - p1p3*2D0))) + 
     &          kp4*(ampi**2*p2p3*(p2p3 + p2p4 + p3p4) - 
     &             ampi**4*p2p4*1D0 + 
     &             p1p3*(ampi**2*p2p4 - 
     &                p2p3*1D0*(ampi**2 + p2p3*2D0)))) - 
     &       kp3*1D0*(ame**2*ampi**2*kp4**2*2D0 + 
     &          p2p4*2D0*
     &           (-(ampi**2*p1p4*p2p4*1D0) + 
     &             p1p3*((ampi**2 + p1p4)*p2p3 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0)) + 
     &             p1p3**2*(p2p4 - p2p3*2D0)) + 
     &          kp2*(-(ampi**2*p1p4*(ampi**2 + p2p4)*1D0) + 
     &             ampi**2*p1p3*(p1p4 + p3p4 - p2p4*1D0) + 
     &             p1p3**2*(ampi**2 + p2p4*2D0 - p3p4*2D0)) + 
     &          kp4*(-(p1p3**2*2D0*
     &                (p2p3 + ame**2*2D0 + kp2*3D0)) + 
     &             p1p3*(-((ame**2 + p1p2)*p3p4*2D0) + 
     &                ampi**2*(-(p2p4*1D0) + ame**2*2D0) + 
     &                p2p3*
     &                 (ampi**2 + ame**2*2D0 - p1p2*2D0 + 
     &                   p1p4*2D0) + ampi**2*kp2*3D0) + 
     &             ampi**2*
     &              (p3p4*x12 - 
     &                (p2p3 + p2p4)*1D0*(ame**2 - p1p2*1D0) + 
     &                p1p4*
     &                 (p2p4 - p2p3*1D0 + ame**2*2D0 + 
     &                   kp2*3D0))))))/
     &   (kp2*kp3*kp4*(ame**2 + p1p2)) - 
     &  (ame**2*D11aex15ex250x34em1em2*ep2*
     &     (kp4*p1p3*(kp2 + p1p4 + p2p3 - kp1*1D0 - 
     &          p1p3*1D0 - p2p4*1D0) + 
     &       kp3*(p1p4*(kp2 + p1p3 + p2p4 - kp1*1D0 - 
     &             p1p4*1D0 - p2p3*1D0) + 
     &          kp4*(ame**2 - p1p2*1D0)*2D0))*4D0)/
     &   (kp3*kp4*(ame**2 + p1p2)) + 
     &  (D33aex15ex250x34em1em2*ep2*
     &     (kp3**2*(ame**2 + p1p2)*p2p4 + 
     &       kp4*p2p3*(kp4*(ame**2 + p1p2) + ame**2*p1p4 + 
     &          ame**2*p2p3 + kp1*p2p3 - ame**2*p1p3*1D0 - 
     &          ame**2*p2p4*1D0 - kp1*p2p4*1D0 - 
     &          kp1*x25*1D0 + kp2*(ame**2 + p1p3 - p1p4*1D0))
     &        - kp3*1D0*(p2p4*
     &           (ame**2*p1p4 + ame**2*p2p3 + kp1*p2p3 + 
     &             kp1*x25 - ame**2*p1p3*1D0 - 
     &             ame**2*p2p4*1D0 - kp1*p2p4*1D0 - 
     &             kp2*1D0*(ame**2 + p1p4 - p1p3*1D0)) + 
     &          kp4*(ame**2*p2p3 + p1p2*p2p3 + 
     &             (ame**2 + p1p2)*p2p4 + ame**4*2D0 + 
     &             ame**2*kp1*2D0 - ame**2*kp2*2D0 - 
     &             ame**2*p1p2*2D0)))*4D0)/
     &   (kp3*kp4*(ame**2 + p1p2)) - 
     &  (D12aex15ex250x34em1em2*ep2*
     &     ((ame**2 + p1p2)*(kp3 - kp4*1D0)*
     &        (kp3*p1p4 - kp4*p1p3*1D0) + 
     &       kp2*(kp3*p1p4**2 - (kp3 + kp4)*p1p3*p1p4*1D0 + 
     &          kp4*(p1p3**2 - ame**2*kp3*2D0)) + 
     &       kp1*(kp4*p1p3*(p2p3 - p2p4*1D0 + kp2*2D0) + 
     &          kp3*(ame**2*kp4*2D0 + 
     &             p1p4*(p2p4 - p2p3*1D0 + kp2*2D0))))*4D0)/
     &   (kp3*kp4*(ame**2 + p1p2)) + 
     &  (D23aex15ex250x34em1em2*ep2*
     &     ((ame**2 + p1p2)*(kp3 - kp4*1D0)*
     &        (kp3*p2p4 - kp4*p2p3*1D0) + 
     &       kp2*(kp3*p2p4*(p1p4 - p1p3*1D0) + 
     &          kp4*(p2p3*(p1p3 - p1p4*1D0) + ame**2*kp3*2D0))
     &         + kp1*(kp4*p2p3*(p2p3 - p2p4*1D0 + kp2*2D0) + 
     &          kp3*(-(ame**2*kp4*2D0) + 
     &             p2p4*(p2p4 - p2p3*1D0 + kp2*2D0))))*4D0)/
     &   (kp3*kp4*(ame**2 + p1p2)) - 
     &  (D00aex15ex250x34em1em2*ep2*
     &     (kp4*(kp2*p1p3 - p1p3*p2p4*1D0 + 
     &          ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p3*(kp1 - p1p4*1D0 + p1p3*2D0)) + 
     &       kp3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*2D0*(p1p2 + ame**2*3D0)))*4D0)/
     &   (kp3*kp4*(ame**2 + p1p2)) + 
     &  (C0appx34m1pm2*ep2*2.5D-1*
     &     (kp1**2*(kp2*(ampi**2*kp3 + kp4*p3p4)*2D0 + 
     &          kp3*(-(ampi**2*p2p3*1D0) + 
     &             p2p4*(ampi**2 + p2p4*2D0)) - 
     &          kp4*1D0*(ampi**2*p2p4 + 
     &             p2p3*(ampi**2 - p2p4*2D0 - p3p4*2D0))) + 
     &       kp2*(kp4*(ampi**2*kp4*(ame**2 + p1p2) - 
     &             ampi**2*kp2*p1p4*1D0 + 
     &             p1p3**2*2D0*(kp2 - p2p4*1D0 + p2p3*2D0) + 
     &             p1p3*(ampi**2*kp2 + p3p4*x12 - 
     &                ampi**2*(ame**2 + p1p2)*2D0 - 
     &                p1p4*p2p3*2D0)) - 
     &          kp3*1D0*(kp2*
     &              (p1p3*(ampi**2 - p1p4*2D0) + 
     &                p1p4*(ampi**2 - p3p4*2D0)) + 
     &             2D0*((ame**2 + p1p2)*(kp4 + p1p4)*p3p4 - 
     &                p1p3*(kp4 + p1p4)*p2p4*1D0 + 
     &                kp4*
     &                 (p1p3*(p1p2 + p2p3 - ame**2*1D0) + 
     &                   ame**2*p2p3*2D0) + 
     &                p1p4*
     &                 (-(ampi**2*(ame**2 + p1p2)*1D0) + 
     &                   p2p3*(-(p1p4*1D0) + p1p3*2D0)))) + 
     &          kp3**2*(ampi**2*(ame**2 + p1p2) - 
     &             p1p4*(ame**2 + p1p2 + p2p4 - p2p3*1D0)*
     &              2D0 - ame**2*kp4*4D0)) + 
     &       kp1*(kp3**2*
     &           (ampi**2*(ame**2 + p1p2) - kp2*p2p4*2D0) + 
     &          kp4*(ampi**2*kp2**2*2D0 + 
     &             kp4*(ampi**2*(ame**2 + p1p2) - 
     &                p2p3*(ame**2 + p1p2 + p1p3 - p1p4*1D0)*
     &                 2D0) + 
     &             p2p3*2D0*
     &              (p1p3*p2p4 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0) + 
     &                p1p4*(p2p3 - p2p4*2D0)) + 
     &             kp2*(ampi**2*(p2p3 - p2p4*1D0) + 
     &                p1p3*(ampi**2 - kp4*2D0 + p2p3*2D0) + 
     &                p1p4*(ampi**2 + p2p3*2D0 - p3p4*2D0)))
     &           + kp3*(ampi**2*kp2*(p1p4 - p1p3*1D0) + 
     &             kp2**2*p3p4*2D0 + 
     &             kp2*(-(p2p3*p3p4*2D0) + 
     &                (p2p3 + p2p4)*(ampi**2 + p1p4*2D0)) - 
     &             p2p4*2D0*
     &              (p1p3*p2p4 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0) + 
     &                p1p4*(p2p3 - p2p4*2D0)) - 
     &             kp4*2D0*
     &              ((ame**2 + p1p2)*p3p4 - 
     &                p2p4*1D0*(ame**2 + p1p3 - p1p2*1D0) + 
     &                kp2*p2p3*3D0 + 
     &                p1p4*(p2p4 + ame**2*2D0 + kp2*3D0)) - 
     &             ame**2*kp4**2*4D0))))/
     &   (kp1*kp2*kp3*kp4*(ame**2 + p1p2)) - 
     &  (D13aex15ex250x34em1em2*ep2*4D0*
     &     (kp3**2*(ame**2 + p1p2)*p1p4 + 
     &       kp4*(kp4*(ame**2 + p1p2)*p1p3 + 
     &          ame**2*p1p3*p1p4 + ame**2*kp1*p2p3 + 
     &          kp1*p1p3*p2p3 + ame**2*p2p3*p2p4 - 
     &          ame**2*p1p3**2*1D0 - ame**2*p1p4*p2p3*1D0 - 
     &          ame**2*p2p3**2*1D0 - ame**2*p1p3*p2p4*1D0 - 
     &          kp1*p1p3*p2p4*1D0 - kp1*p1p3*x25*1D0 + 
     &          kp2*(-(ame**2*p2p3*1D0) + 
     &             p1p3*(ame**2 + p1p3 - p1p4*1D0)) + 
     &          ame**2*p1p3*p2p3*2D0) + 
     &       kp3*(ame**2*p1p3*p1p4 + ame**2*kp1*p2p4 + 
     &          kp1*p1p4*p2p4 + ame**2*p2p3*p2p4 - 
     &          ame**2*p1p4**2*1D0 - ame**2*p1p4*p2p3*1D0 - 
     &          kp1*p1p4*p2p3*1D0 - ame**2*p1p3*p2p4*1D0 - 
     &          ame**2*p2p4**2*1D0 - kp1*p1p4*x25*1D0 + 
     &          kp2*(-(ame**2*p2p4*1D0) + 
     &             p1p4*(ame**2 + p1p4 - p1p3*1D0)) + 
     &          ame**2*p1p4*p2p4*2D0 + 
     &          kp4*(-(ame**2*p1p3*1D0) - p1p2*p1p3*1D0 - 
     &             (ame**2 + p1p2)*p1p4*1D0 + 
     &             ame**2*kp1*2D0 - ame**2*kp2*2D0 + 
     &             ame**4*4D0 - ame**2*p1p2*4D0))))/
     &   (kp3*kp4*(ame**2 + p1p2)) - 
     &  (C0aex34x15em2m1*ep2*5.D-1*
     &     (kp3**2*2D0*(ampi**2*(ame**2 + p1p2) - 
     &          p1p4*1D0*(ame**2 + p1p2 + p2p4 - p2p3*1D0) - 
     &          ame**2*kp4*2D0) - 
     &       kp3*1D0*(-(p1p4*p3p4*x12*1D0) + 
     &          2D0*(ampi**2*kp2*p1p3 + ame**2*kp4**2*2D0 + 
     &             p1p4*((ame**2 + p1p2)*p3p4 - 
     &                kp2*(p1p3 + p1p4 + p3p4)*1D0 + 
     &                p1p3*p2p3*2D0 - p1p4*p2p4*2D0) + 
     &             kp4*(p1p4*
     &                 (p1p2 + p2p4 - ame**2*1D0 - p2p3*1D0)
     &                 + p1p3*
     &                 (p1p2 + p2p3 - ame**2*1D0 - p2p4*1D0)
     &                 + ame**2*(p2p3 + p2p4)*2D0 + 
     &                (ame**2 + p1p2)*p3p4*2D0))) + 
     &       kp1*2D0*(kp2*(kp3 + kp4)*(ampi**2 + p3p4) - 
     &          kp4**2*p2p3*1D0 + 
     &          kp3*(p1p4*(p2p3 + p2p4) - p2p3*p3p4*1D0 + 
     &             p2p4*(ampi**2 - kp3*1D0)) + 
     &          kp4*(p2p3*(ampi**2 + p1p3 - kp3*3D0) + 
     &             p2p4*(p1p3 - p3p4*1D0 - kp3*3D0))) + 
     &       kp4*(kp2*(p1p3*(p1p3 + p1p4 + p3p4) - 
     &             ampi**2*p1p4*1D0)*2D0 + 
     &          kp4*(ame**2*ampi**2 + 
     &             p1p2*(ampi**2 - p1p3*1D0) - 
     &             p1p3*1D0*(ame**2 + p2p3 - p2p4*1D0))*2D0 + 
     &          p1p3*(p3p4*(x12 - (ame**2 + p1p2)*2D0) + 
     &             p1p3*p2p3*4D0 - p1p4*p2p4*4D0))))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) - 
     &  (C0aex34x25em1m2*ep2*5.D-1*
     &     (kp1*(kp2*(kp3 + kp4)*(ampi**2 + p3p4) + 
     &          kp4*(p2p3*(p2p3 + p2p4 + p3p4) - 
     &             ampi**2*p2p4*1D0) + 
     &          kp3*(p2p4*(p2p4 + p3p4) + 
     &             p2p3*(p2p4 - ampi**2*1D0)))*2D0 + 
     &       kp3**2*2D0*(ampi**2*(ame**2 + p1p2) - 
     &          kp2*p1p4*1D0 - 
     &          p2p4*1D0*(ame**2 + p1p2 + p1p4 - p1p3*1D0) - 
     &          ame**2*kp4*2D0) - 
     &       kp3*2D0*(kp2*p1p3*p3p4 - ampi**2*kp2*p1p4*1D0 - 
     &          kp2*p1p3*p2p4*1D0 - kp2*p1p4*p2p4*1D0 + 
     &          ame**2*kp4**2*2D0 + p1p3*p2p3*p2p4*2D0 - 
     &          p1p4*p2p4**2*2D0 + 
     &          kp4*(-((p2p3 + p2p4)*1D0*
     &                (ame**2 - p1p2*1D0)) + 
     &             (ame**2 + p1p2)*p3p4*2D0 + 
     &             p1p4*(p2p4 - p2p3*1D0 + ame**2*2D0 + 
     &                kp2*3D0) + 
     &             p1p3*(p2p3 - p2p4*1D0 + ame**2*2D0 + 
     &                kp2*3D0))) + 
     &       kp4*(kp4*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp2*p1p3*1D0 - 
     &             p2p3*1D0*(ame**2 + p1p2 + p1p3 - p1p4*1D0))
     &            *2D0 + 
     &          kp2*(p1p3*(ampi**2 + p2p3) + 
     &             p1p4*(p2p3 - p3p4*1D0))*2D0 + 
     &          p2p3*(p3p4*(x12 - (ame**2 + p1p2)*2D0) + 
     &             p1p3*p2p3*4D0 - p1p4*p2p4*4D0))))/
     &   (kp2*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D1aex12px23ex45em1m2p*ep3*2D0*
     &     (kp1**2*p2p3*(ampi**2*kp3 + kp4*p3p4) + 
     &       ampi**2*kp3*
     &        (p1p4**2*p2p3 + 
     &          p2p3*((ame**2 + p1p2)*(kp3 + p3p4) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 - 
     &             p1p3*(kp2 + p2p4)*1D0) + 
     &          p1p4*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             (ame**2 + p1p2)*(kp3 + p3p4)*1D0 - 
     &             p2p3*1D0*(p2p3 + p1p3*2D0) + 
     &             kp2*(p1p3 + p2p3*2D0) + 
     &             p2p4*(p1p3 + p2p3*2D0))) + 
     &       kp4**2*2D0*(p1p3**2*
     &           (p2p4 - kp2*1D0 - p2p3*2D0) + 
     &          p1p3*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             (ame**2 + p1p2)*p3p4*1D0 + 
     &             kp3*(p1p2 - ame**2*1D0) + 
     &             p2p3*(p1p4 + p2p4 - kp2*1D0 - p2p3*2D0)) + 
     &          p2p3*(ampi**2*(ame**2 + p1p2) + p1p4*p2p3 - 
     &             (ame**2 + p1p2)*p3p4*1D0 - 
     &             kp3*1D0*(p1p2 + ame**2*3D0))) + 
     &       kp4*((ame**2 + p1p2)*p3p4**2*(p1p4 - p2p3*1D0) + 
     &          kp3*2D0*(p2p3*
     &              (ampi**2*(ame**2 + p1p2) + p1p3*p2p4) - 
     &             p1p4**2*p2p3*2D0 + 
     &             p1p4*(kp3*(ame**2 + p1p2) + p2p3**2 - 
     &                ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &                p1p3*p2p3*2D0 - 
     &                kp2*1D0*(p1p3 + p2p3*2D0) - 
     &                p2p4*2D0*(p1p3 + p2p3*2D0))) + 
     &          p3p4*(-(p1p4**2*p2p3*1D0) + 
     &             p2p3*(ame**2*ampi**2 + ampi**2*p1p2 + 
     &                p1p3*(kp2 + p2p4) - 
     &                kp3*(ame**2 + p1p2)*3D0) + 
     &             p1p4*(kp2*p1p3 + p2p3**2 - 
     &                ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &                p1p3*p2p3*2D0 - 
     &                p2p4*1D0*(p1p3 + p2p3*2D0) + 
     &                kp3*(ame**2 + p1p2)*3D0))) + 
     &       kp1*(-(ampi**2*kp3**2*(ame**2 + p1p2)*1D0) + 
     &          kp4*(-(kp4*p2p3*(p1p3 + p2p3)*2D0) + 
     &             p3p4*(p2p3**2 + 
     &                ame**2*(p3p4 - ampi**2*1D0) + 
     &                p1p2*(p3p4 - ampi**2*1D0) + 
     &                kp2*(p1p3 + p2p3*2D0) + 
     &                p1p3*(-(p2p4*1D0) + p2p3*2D0))) + 
     &          kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &             ampi**2*kp2*p1p3 + ampi**2*p1p3*p2p4 - 
     &             ampi**2*p2p3**2*1D0 - 
     &             (ame**2 + p1p2)*p3p4*1D0*
     &              (ampi**2 - kp4*1D0) + 
     &             ampi**2*kp2*p2p3*2D0 - 
     &             ampi**2*p1p3*p2p3*2D0 + 
     &             ampi**2*p1p4*p2p3*2D0 - 
     &             kp4*p1p3*p2p4*2D0 + 
     &             ampi**2*p2p3*p2p4*2D0 - 
     &             kp4*p1p4*p2p3*4D0 - kp4*p2p3*p2p4*4D0))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)) + 
     &  (C1apx34ppm2m1*ep2*2.5D-1*
     &     (-(kp1**2*(kp2*(kp3 + kp4)*(ampi**2 + p3p4) + 
     &            kp4*(p2p3*(p2p3 + p2p4 + p3p4) - 
     &               ampi**2*p2p4*1D0) + 
     &            kp3*(p2p4*(p2p4 + p3p4) + 
     &               p2p3*(p2p4 - ampi**2*1D0)))*2D0) - 
     &       kp2*1D0*(-(kp3**2*2D0*
     &             (ampi**2*(ame**2 + p1p2) - 
     &               p1p4*1D0*
     &                (ame**2 + p1p2 + p2p4 - p2p3*1D0) - 
     &               ame**2*kp4*2D0)) - 
     &          kp4*1D0*(kp2*
     &              (p1p3*(p1p3 + p1p4 + p3p4) - 
     &                ampi**2*p1p4*1D0)*2D0 + 
     &             kp4*(ame**2*ampi**2 + 
     &                p1p2*(ampi**2 - p1p3*1D0) - 
     &                p1p3*1D0*(ame**2 + p2p3 - p2p4*1D0))*2D0
     &               + p1p3*
     &              (p3p4*(x12 - (ame**2 + p1p2)*2D0) + 
     &                p1p3*p2p3*4D0 - p1p4*p2p4*4D0)) + 
     &          kp3*(kp2*
     &              (-(p1p4*(p1p4 + p3p4)*1D0) + 
     &                p1p3*(ampi**2 - p1p4*1D0))*2D0 + 
     &             kp4*2D0*
     &              (p1p4*
     &                 (p1p2 + p2p4 - ame**2*1D0 - p2p3*1D0)
     &                 + p1p3*
     &                 (p1p2 + p2p3 - ame**2*1D0 - p2p4*1D0)
     &                 + ame**2*(p2p3 + p2p4)*2D0 + 
     &                (ame**2 + p1p2)*p3p4*2D0) + 
     &             ame**2*kp4**2*4D0 + 
     &             p1p4*(p3p4*
     &                 (-(x12*1D0) + (ame**2 + p1p2)*2D0) + 
     &                p1p3*p2p3*4D0 - p1p4*p2p4*4D0))) + 
     &       kp1*(kp2**2*(kp3 + kp4)*(ampi**2 + p3p4)*2D0 + 
     &          kp3**2*2D0*
     &           (-(ampi**2*(ame**2 + p1p2)*1D0) + 
     &             p2p4*(ame**2 + p1p2 + p1p4 - p1p3*1D0) + 
     &             ame**2*kp4*2D0) + 
     &          kp3*2D0*(ame**2*kp4**2*2D0 + 
     &             p2p4*(p1p3*p2p3 - p1p4*p2p4*1D0)*2D0 + 
     &             kp4*(-((p2p3 + p2p4)*1D0*
     &                   (ame**2 - p1p2*1D0)) + 
     &                (ame**2 + p1p2)*p3p4*2D0 + 
     &                p1p4*(p2p4 - p2p3*1D0 + ame**2*2D0) + 
     &                p1p3*(p2p3 - p2p4*1D0 + ame**2*2D0))) - 
     &          kp2*2D0*(kp4**2*(p2p3 - p1p3*1D0) + 
     &             kp3*(p2p4*(kp3 + p1p3 - ampi**2*1D0) + 
     &                p3p4*(p2p3 - p1p3*1D0) + 
     &                p1p4*(ampi**2 - kp3*1D0 - p2p3*1D0)) + 
     &             kp4*(p2p3*(p1p4 - ampi**2*1D0) + 
     &                p3p4*(p2p4 - p1p4*1D0) + 
     &                kp3*(p2p3 + p2p4 - p1p4*1D0)*3D0 + 
     &                p1p3*(ampi**2 - p2p4*1D0 - kp3*3D0))) + 
     &          kp4*(kp4*
     &              (-(ampi**2*(ame**2 + p1p2)*1D0) + 
     &                p2p3*(ame**2 + p1p2 + p1p3 - p1p4*1D0))*
     &              2D0 + 
     &             p2p3*(p3p4*
     &                 (-(x12*1D0) + (ame**2 + p1p2)*2D0) - 
     &                p1p3*p2p3*4D0 + p1p4*p2p4*4D0)))))/
     &   (kp1*kp2*kp3*kp4*(ame**2 + p1p2)) + 
     &  (C1apx34ppm1m2*ep2*2.5D-1*
     &     (kp1**2*(kp2*(kp3 + kp4)*(ampi**2 + p3p4) + 
     &          kp4*(p2p3*(p2p3 + p2p4 + p3p4) - 
     &             ampi**2*p2p4*1D0) + 
     &          kp3*(p2p4*(p2p4 + p3p4) + 
     &             p2p3*(p2p4 - ampi**2*1D0)))*2D0 + 
     &       kp2*(-(kp3**2*2D0*
     &             (ampi**2*(ame**2 + p1p2) - 
     &               p1p4*1D0*
     &                (ame**2 + p1p2 + p2p4 - p2p3*1D0) - 
     &               ame**2*kp4*2D0)) - 
     &          kp4*1D0*(kp2*
     &              (p1p3*(p1p3 + p1p4 + p3p4) - 
     &                ampi**2*p1p4*1D0)*2D0 + 
     &             kp4*(ame**2*ampi**2 + 
     &                p1p2*(ampi**2 - p1p3*1D0) - 
     &                p1p3*1D0*(ame**2 + p2p3 - p2p4*1D0))*2D0
     &               + p1p3*
     &              (p3p4*(x12 - (ame**2 + p1p2)*2D0) + 
     &                p1p3*p2p3*4D0 - p1p4*p2p4*4D0)) + 
     &          kp3*(kp2*
     &              (-(p1p4*(p1p4 + p3p4)*1D0) + 
     &                p1p3*(ampi**2 - p1p4*1D0))*2D0 + 
     &             kp4*2D0*
     &              (p1p4*
     &                 (p1p2 + p2p4 - ame**2*1D0 - p2p3*1D0)
     &                 + p1p3*
     &                 (p1p2 + p2p3 - ame**2*1D0 - p2p4*1D0)
     &                 + ame**2*(p2p3 + p2p4)*2D0 + 
     &                (ame**2 + p1p2)*p3p4*2D0) + 
     &             ame**2*kp4**2*4D0 + 
     &             p1p4*(p3p4*
     &                 (-(x12*1D0) + (ame**2 + p1p2)*2D0) + 
     &                p1p3*p2p3*4D0 - p1p4*p2p4*4D0))) - 
     &       kp1*1D0*(kp2**2*(kp3 + kp4)*(ampi**2 + p3p4)*
     &           2D0 + kp3**2*2D0*
     &           (-(ampi**2*(ame**2 + p1p2)*1D0) + 
     &             p2p4*(ame**2 + p1p2 + p1p4 - p1p3*1D0) + 
     &             ame**2*kp4*2D0) + 
     &          kp3*2D0*(ame**2*kp4**2*2D0 + 
     &             p2p4*(p1p3*p2p3 - p1p4*p2p4*1D0)*2D0 + 
     &             kp4*(-((p2p3 + p2p4)*1D0*
     &                   (ame**2 - p1p2*1D0)) + 
     &                (ame**2 + p1p2)*p3p4*2D0 + 
     &                p1p4*(p2p4 - p2p3*1D0 + ame**2*2D0) + 
     &                p1p3*(p2p3 - p2p4*1D0 + ame**2*2D0))) - 
     &          kp2*2D0*(kp4**2*(p2p3 - p1p3*1D0) + 
     &             kp3*(p2p4*(kp3 + p1p3 - ampi**2*1D0) + 
     &                p3p4*(p2p3 - p1p3*1D0) + 
     &                p1p4*(ampi**2 - kp3*1D0 - p2p3*1D0)) + 
     &             kp4*(p2p3*(p1p4 - ampi**2*1D0) + 
     &                p3p4*(p2p4 - p1p4*1D0) + 
     &                kp3*(p2p3 + p2p4 - p1p4*1D0)*3D0 + 
     &                p1p3*(ampi**2 - p2p4*1D0 - kp3*3D0))) + 
     &          kp4*(kp4*
     &              (-(ampi**2*(ame**2 + p1p2)*1D0) + 
     &                p2p3*(ame**2 + p1p2 + p1p3 - p1p4*1D0))*
     &              2D0 + 
     &             p2p3*(p3p4*
     &                 (-(x12*1D0) + (ame**2 + p1p2)*2D0) - 
     &                p1p3*p2p3*4D0 + p1p4*p2p4*4D0)))))/
     &   (kp1*kp2*kp3*kp4*(ame**2 + p1p2)) - 
     &  (D1aex12px24ex35em1m2p*ep3*2D0*
     &     (kp1**2*p2p4*(ampi**2*kp4 + kp3*p3p4) - 
     &       kp4**2*(ame**2 + p1p2)*1D0*
     &        (-(ampi**2*p2p4*1D0) + p1p3*(ampi**2 - kp3*2D0))
     &         + kp3*(kp3*(p1p4 + p2p4)*2D0*
     &           (p1p3*p2p4 + ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p1p4*(p2p3 - kp2*1D0 - p2p4*2D0)) + 
     &          p3p4*(-(p1p3**2*p2p4*1D0) + 
     &             p2p4*(p1p4*(kp2 + p2p3) + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0)) + 
     &             p1p3*(p2p4**2 + 
     &                ame**2*(p3p4 - ampi**2*1D0) + 
     &                p1p2*(p3p4 - ampi**2*1D0) + 
     &                p1p4*(kp2 - p2p3*1D0) + 
     &                p2p4*(p1p4 - p2p3*1D0)*2D0))) + 
     &       kp1*(-(ampi**2*kp4**2*(ame**2 + p1p2)*1D0) + 
     &          kp3*(-(kp3*p2p4*(p1p4 + p2p4)*2D0) + 
     &             p3p4*(p2p4**2 + 
     &                ame**2*(p3p4 - ampi**2*1D0) + 
     &                p1p2*(p3p4 - ampi**2*1D0) + 
     &                p1p4*(kp2 - p2p3*1D0) + 
     &                (kp2 + p1p4)*p2p4*2D0)) + 
     &          kp4*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &             ampi**2*kp2*p1p4 + ampi**2*p1p4*p2p3 - 
     &             ampi**2*p2p4**2*1D0 - 
     &             (ame**2 + p1p2)*p3p4*1D0*
     &              (ampi**2 - kp3*1D0) - kp3*p1p4*p2p3*2D0 + 
     &             ampi**2*kp2*p2p4*2D0 + 
     &             ampi**2*p1p3*p2p4*2D0 - 
     &             ampi**2*p1p4*p2p4*2D0 + 
     &             ampi**2*p2p3*p2p4*2D0 - 
     &             kp3*p1p3*p2p4*4D0 - kp3*p2p3*p2p4*4D0)) + 
     &       kp4*(-(ampi**2*p2p4*1D0*
     &             (p1p4*(kp2 + p2p3) + 
     &               ame**2*(ampi**2 - p3p4*1D0) + 
     &               p1p2*(ampi**2 - p3p4*1D0))) + 
     &          kp3*p2p4*
     &           (ampi**2*(ame**2 + p1p2)*2D0 + 
     &             p1p4*p2p3*2D0 - (ame**2 + p1p2)*p3p4*3D0)
     &           - kp3**2*2D0*
     &           (p1p4*(ame**2 - p1p2*1D0) + 
     &             p2p4*(p1p2 + ame**2*3D0)) + 
     &          p1p3**2*p2p4*(ampi**2 - kp3*4D0) + 
     &          p1p3*(ampi**2*
     &              (p1p4*(kp2 + p2p3) - p2p4**2*1D0 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0) + 
     &                p2p4*(kp2 + p2p3 - p1p4*1D0)*2D0) + 
     &             kp3*(p2p4**2*2D0 - 
     &                2D0*
     &                 (ame**2*ampi**2 + ampi**2*p1p2 + 
     &                   p1p4*(kp2 + p2p3*2D0)) + 
     &                (ame**2 + p1p2)*p3p4*3D0 - 
     &                p2p4*(kp2 - p1p4*1D0 + p2p3*2D0)*4D0))))
     &     )/(kp3**2*kp4*(ame**2 + p1p2)) + 
     &  (C0ax12px35m1m2p*ep3*
     &     (ampi**2*kp4*(-(1D0*(kp1 + p1p3 - p1p4*1D0)*
     &             (kp2 + p2p3 - p2p4*1D0)) + 
     &          (ame**2 + p1p2)*(ampi**2 - kp4*1D0 - p3p4*1D0)
     &          ) + kp3*((ame**2 + p1p2)*p3p4**2 - 
     &          p3p4*1D0*
     &           (ame**2*ampi**2 + kp1*kp2 + ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0) - 
     &             kp4*(ame**2 + p1p2)*2D0) - 
     &          kp4*1D0*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp2*p1p3*2D0 + p1p3*p2p4*2D0 - 
     &             p2p3*2D0*(kp1 - p1p4*1D0 + p1p3*2D0))) + 
     &       kp3**2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*(p1p2 + ame**2*5D0))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)) - 
     &  (D0aex12x35x14epem1m2p*ep3*p1p4*4D0*
     &     (ampi**2*kp4*(-(1D0*(kp1 + p1p3 - p1p4*1D0)*
     &             (kp2 + p2p3 - p2p4*1D0)) + 
     &          (ame**2 + p1p2)*(ampi**2 - kp4*1D0 - p3p4*1D0)
     &          ) + kp3*((ame**2 + p1p2)*p3p4**2 - 
     &          p3p4*1D0*
     &           (ame**2*ampi**2 + kp1*kp2 + ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0) - 
     &             kp4*(ame**2 + p1p2)*2D0) - 
     &          kp4*1D0*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp2*p1p3*2D0 + p1p3*p2p4*2D0 - 
     &             p2p3*2D0*(kp1 - p1p4*1D0 + p1p3*2D0))) + 
     &       kp3**2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*(p1p2 + ame**2*5D0))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)) + 
     &  (D0aex12px24ex35em1m2p*ep3*p2p4*4D0*
     &     (ampi**2*kp4*(-(1D0*(kp1 + p1p3 - p1p4*1D0)*
     &             (kp2 + p2p3 - p2p4*1D0)) + 
     &          (ame**2 + p1p2)*(ampi**2 - kp4*1D0 - p3p4*1D0)
     &          ) + kp3*((ame**2 + p1p2)*p3p4**2 - 
     &          p3p4*1D0*
     &           (ame**2*ampi**2 + kp1*kp2 + ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0) - 
     &             kp4*(ame**2 + p1p2)*2D0) - 
     &          kp4*1D0*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp2*p1p3*2D0 + p1p3*p2p4*2D0 - 
     &             p2p3*2D0*(kp1 - p1p4*1D0 + p1p3*2D0))) + 
     &       kp3**2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*(p1p2 + ame**2*5D0))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)) + 
     &  (D3aex12px24ex35em1m2p*ep3*2D0*
     &     (-(ampi**2*kp4*1D0*
     &          (ampi**2*p1p3*p2p3 + ame**2*kp4*p2p4 + 
     &            kp4*p1p2*p2p4 + (ame**2 + p1p2)*p2p4*p3p4 - 
     &            ame**2*ampi**2*p2p4*1D0 - 
     &            ampi**2*p1p2*p2p4*1D0 - 
     &            ampi**2*p1p4*p2p4*1D0 - 
     &            p1p4*p2p3*p2p4*1D0 - p1p3*p2p4**2*1D0 + 
     &            kp2*(ampi**2*p1p3 - p1p4*p2p4*1D0) + 
     &            kp1*(ampi**2*(kp2 + p2p3) - p2p4**2*1D0) + 
     &            p1p4*p2p4**2*2D0)) + 
     &       kp3*((ame**2 + p1p2)*p2p4*p3p4**2 + 
     &          p3p4*(ampi**2*
     &              (kp4*(ame**2 + p1p2) + p1p3*p2p3 - 
     &                kp2*(kp1 + p1p4)*1D0) - 
     &             p2p4*1D0*
     &              (ame**2*ampi**2 + ampi**2*kp1 + 
     &                ampi**2*p1p2 - 
     &                kp4*(ame**2 + p1p2)*1D0 + 
     &                p1p4*(ampi**2 + p2p3 - kp2*1D0)) + 
     &             p2p4**2*(kp1 - p1p3*1D0 + p1p4*2D0)) - 
     &          kp4*1D0*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &             p1p4*p2p3*p2p4*2D0 + p1p3*p2p4**2*2D0 - 
     &             ampi**2*2D0*
     &              (kp2*p1p3 + p2p3*(kp1 + p1p3*2D0)))) + 
     &       ampi**2*kp3**2*
     &        (ame**2*(p3p4 - ampi**2*1D0) + 
     &          p1p2*(p3p4 - ampi**2*1D0) + 
     &          p2p4*(kp1 - p1p3*1D0) + 
     &          p1p4*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*(p1p2 + ame**2*5D0))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)) + 
     &  (C0ax12px35m2m1p*ep3*
     &     (ampi**2*kp4*((kp1 + p1p3 - p1p4*1D0)*
     &           (kp2 + p2p3 - p2p4*1D0) - 
     &          (ame**2 + p1p2)*1D0*
     &           (ampi**2 - kp4*1D0 - p3p4*1D0)) + 
     &       kp3*(-((ame**2 + p1p2)*p3p4**2*1D0) + 
     &          p3p4*(ame**2*ampi**2 + kp1*kp2 + 
     &             ampi**2*p1p2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0) - 
     &             kp4*(ame**2 + p1p2)*2D0) + 
     &          kp4*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &             kp2*p1p3*2D0 + p1p3*p2p4*2D0 - 
     &             p2p3*2D0*(kp1 - p1p4*1D0 + p1p3*2D0))) + 
     &       kp3**2*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &          (ame**2 + p1p2)*p3p4*1D0 + 
     &          p2p4*(p1p3 - kp1*1D0) + 
     &          p1p4*(p2p3 - kp2*1D0 - p2p4*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*5D0))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)) + 
     &  (D3aex12x35x14epem1m2p*ep3*2D0*
     &     (ampi**2*kp4*(ame**2*kp4*p1p4 + kp4*p1p2*p1p4 + 
     &          ampi**2*p1p3*p2p3 + 
     &          (ame**2 + p1p2)*p1p4*p3p4 - 
     &          ame**2*ampi**2*p1p4*1D0 - 
     &          ampi**2*p1p2*p1p4*1D0 - p1p4**2*p2p3*1D0 - 
     &          ampi**2*p1p4*p2p4*1D0 - p1p3*p1p4*p2p4*1D0 + 
     &          kp2*(ampi**2*p1p3 - p1p4**2*1D0) + 
     &          kp1*(ampi**2*(kp2 + p2p3) - p1p4*p2p4*1D0) + 
     &          p1p4**2*p2p4*2D0) + 
     &       kp3*(p3p4*(-(ampi**2*p1p3*p2p3*1D0) + 
     &             kp1*(ampi**2*kp2 + 
     &                p2p4*(ampi**2 - p1p4*1D0)) + 
     &             p1p4*(ampi**2*kp2 + 
     &                (ampi**2 + p1p3)*p2p4 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0)) + 
     &             p1p4**2*(p2p3 - kp2*1D0 - p2p4*2D0)) + 
     &          kp4*(ame**2*ampi**4 - 
     &             ame**2*(ampi**2 + p1p4)*p3p4*1D0 + 
     &             p1p2*(ampi**4 - 
     &                (ampi**2 + p1p4)*p3p4*1D0) - 
     &             ampi**2*kp2*p1p3*2D0 + 
     &             p1p3*p1p4*p2p4*2D0 + 
     &             p2p3*2D0*
     &              (p1p4**2 - ampi**2*1D0*(kp1 + p1p3*2D0))))
     &         + ampi**2*kp3**2*
     &        (ame**2*ampi**2 + ampi**2*p1p2 - 
     &          (ame**2 + p1p2)*p3p4*1D0 + 
     &          p2p4*(p1p3 - kp1*1D0) + 
     &          p1p4*(p2p3 - kp2*1D0 - p2p4*2D0) - 
     &          kp4*1D0*(p1p2 + ame**2*5D0))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)) + 
     &  (C0ax12px45m1m2p*ep3*
     &     (ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*(p3p4*(kp1*kp2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0) + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0)) + 
     &          kp4*(p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0))) + 
     &       kp3*(ampi**2*
     &           (ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) - 
     &             1D0*(p1p3 - p1p4*1D0)*
     &              (kp2 + p2p4 - p2p3*1D0)) + 
     &          kp1*(ampi**2*(kp2 - p2p3*1D0) + 
     &             p2p4*(ampi**2 - kp4*2D0)) + 
     &          kp4*(p1p3*p2p4*2D0 - 
     &             p1p4*2D0*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &             ame**2*(ampi**2 - p3p4*2D0) + 
     &             p1p2*(ampi**2 - p3p4*2D0)) - 
     &          kp4**2*1D0*(p1p2 + ame**2*5D0))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)) - 
     &  (C0ax12px45m2m1p*ep3*1D0*
     &     (ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*(p3p4*(kp1*kp2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0) + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0)) + 
     &          kp4*(p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0))) + 
     &       kp3*(ampi**2*
     &           (ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) - 
     &             1D0*(p1p3 - p1p4*1D0)*
     &              (kp2 + p2p4 - p2p3*1D0)) + 
     &          kp1*(ampi**2*(kp2 - p2p3*1D0) + 
     &             p2p4*(ampi**2 - kp4*2D0)) + 
     &          kp4*(p1p3*p2p4*2D0 - 
     &             p1p4*2D0*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &             ame**2*(ampi**2 - p3p4*2D0) + 
     &             p1p2*(ampi**2 - p3p4*2D0)) - 
     &          kp4**2*1D0*(p1p2 + ame**2*5D0))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)) - 
     &  (D0aex12x45x13epem1m2p*ep3*p1p3*4D0*
     &     (ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*(p3p4*(kp1*kp2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0) + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0)) + 
     &          kp4*(p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0))) + 
     &       kp3*(ampi**2*
     &           (ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) - 
     &             1D0*(p1p3 - p1p4*1D0)*
     &              (kp2 + p2p4 - p2p3*1D0)) + 
     &          kp1*(ampi**2*(kp2 - p2p3*1D0) + 
     &             p2p4*(ampi**2 - kp4*2D0)) + 
     &          kp4*(p1p3*p2p4*2D0 - 
     &             p1p4*2D0*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &             ame**2*(ampi**2 - p3p4*2D0) + 
     &             p1p2*(ampi**2 - p3p4*2D0)) - 
     &          kp4**2*1D0*(p1p2 + ame**2*5D0))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)) + 
     &  (D0aex12px23ex45em1m2p*ep3*p2p3*4D0*
     &     (ampi**2*kp3**2*(ame**2 + p1p2) + 
     &       kp4*(p3p4*(kp1*kp2 - 
     &             1D0*(p1p3 - p1p4*1D0)*(p2p3 - p2p4*1D0) + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0)) + 
     &          kp4*(p1p3*p2p4 - kp2*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0))) + 
     &       kp3*(ampi**2*
     &           (ame**2*(p3p4 - ampi**2*1D0) + 
     &             p1p2*(p3p4 - ampi**2*1D0) - 
     &             1D0*(p1p3 - p1p4*1D0)*
     &              (kp2 + p2p4 - p2p3*1D0)) + 
     &          kp1*(ampi**2*(kp2 - p2p3*1D0) + 
     &             p2p4*(ampi**2 - kp4*2D0)) + 
     &          kp4*(p1p3*p2p4*2D0 - 
     &             p1p4*2D0*(kp2 - p2p3*1D0 + p2p4*2D0) + 
     &             ame**2*(ampi**2 - p3p4*2D0) + 
     &             p1p2*(ampi**2 - p3p4*2D0)) - 
     &          kp4**2*1D0*(p1p2 + ame**2*5D0))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)) - 
     &  (D3aex12x45x13epem1m2p*ep3*2D0*
     &     (ampi**2*kp3**2*(ame**2 + p1p2)*p1p3 + 
     &       kp4*(-((ame**2 + p1p2)*p1p3*p3p4**2*1D0) + 
     &          ampi**2*kp4*
     &           (ame**2*ampi**2 + ampi**2*p1p2 + p1p4*p2p3 - 
     &             kp2*p1p3*1D0 + p1p3*(p2p4 - p2p3*2D0)) + 
     &          p3p4*(p1p3*
     &              (ampi**2*(ame**2 + kp2 + p1p2) + 
     &                (ampi**2 + p1p4)*p2p3) - 
     &             ampi**2*kp4*(ame**2 + p1p2)*1D0 - 
     &             ampi**2*p1p4*p2p4*1D0 + 
     &             p1p3**2*(p2p4 - kp2*1D0 - p2p3*2D0))) + 
     &       kp1*(kp4*p3p4*
     &           (ampi**2*kp2 + p2p3*(ampi**2 - p1p3*1D0)) + 
     &          ampi**2*(-(kp4**2*p2p3*1D0) + 
     &             kp3*(ampi**2*kp2 - p1p3*p2p3*1D0 + 
     &                p2p4*(ampi**2 - kp4*2D0)))) + 
     &       kp3*(ampi**2*
     &           (ampi**2*p1p4*(kp2 + p2p4) + 
     &             p1p3*((ame**2 + p1p2)*p3p4 - 
     &                ampi**2*(ame**2 + p1p2)*1D0 - 
     &                (ampi**2 + p1p4)*p2p3*1D0) - 
     &             p1p3**2*1D0*(kp2 + p2p4 - p2p3*2D0)) + 
     &          kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &             (ame**2 + p1p2)*(ampi**2 + p1p3)*p3p4*
     &              1D0 - ampi**2*kp2*p1p4*2D0 + 
     &             p1p3*p1p4*p2p3*2D0 + p1p3**2*p2p4*2D0 - 
     &             ampi**2*p1p4*p2p4*4D0) - 
     &          ampi**2*kp4**2*1D0*(p1p2 + ame**2*5D0))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)) + 
     &  (D3aex12px23ex45em1m2p*ep3*2D0*
     &     (ampi**2*kp3**2*(ame**2 + p1p2)*p2p3 + 
     &       kp4*(p3p4*(-(ampi**2*p1p4*p2p4*1D0) + 
     &             kp2*p1p3*(ampi**2 - p2p3*1D0) + 
     &             p2p3*(p1p3*(ampi**2 + p2p4) + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0)) + 
     &             p2p3**2*(p1p4 - p1p3*2D0)) + 
     &          ampi**2*kp4*
     &           (p1p4*p2p3 - kp2*p1p3*1D0 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p1p3*(p2p4 - p2p3*2D0))) + 
     &       kp1*(kp4*p3p4*
     &           (ampi**2*kp2 + p2p3*(ampi**2 - p2p3*1D0)) + 
     &          ampi**2*(-(kp4**2*p2p3*1D0) + 
     &             kp3*(ampi**2*kp2 - p2p3**2*1D0 + 
     &                p2p4*(ampi**2 - kp4*2D0)))) + 
     &       kp3*(ampi**2*
     &           (ampi**2*p1p4*p2p4 + 
     &             kp2*(ampi**2*p1p4 - p1p3*p2p3*1D0) + 
     &             p2p3*((ame**2 + p1p2)*p3p4 - 
     &                ampi**2*(ame**2 + p1p2 + p1p3)*1D0 - 
     &                p1p3*p2p4*1D0) + 
     &             p2p3**2*(-(p1p4*1D0) + p1p3*2D0)) + 
     &          kp4*(ame**2*ampi**4 + ampi**4*p1p2 - 
     &             (ame**2 + p1p2)*(ampi**2 + p2p3)*p3p4*
     &              1D0 - ampi**2*kp2*p1p4*2D0 + 
     &             p1p4*p2p3**2*2D0 + p1p3*p2p3*p2p4*2D0 - 
     &             ampi**2*p1p4*p2p4*4D0) - 
     &          ampi**2*kp4**2*1D0*(p1p2 + ame**2*5D0))))/
     &   (kp3*kp4**2*(ame**2 + p1p2)) + 
     &  (D3ax15x34px23epem1m2p*ep2*
     &     (kp3**2*(ame**2*ampi**4 + ampi**2*p1p4*p2p3 - 
     &          ame**2*ampi**2*p1p4*1D0 - 
     &          ame**2*ampi**2*p2p3*1D0 + 
     &          ame**2*p1p4*p2p3*2D0 - p1p4*p2p3**2*2D0 - 
     &          p1p4*p2p4*1D0*(ampi**2 - p2p3*2D0) - 
     &          ame**2*kp4*2D0*(ampi**2 - p2p3*2D0) + 
     &          p1p2*(ampi**4 - ampi**2*p1p4*1D0 - 
     &             p2p3*1D0*(ampi**2 - p1p4*2D0))) + 
     &       kp4*(ampi**2*kp4*
     &           (ame**2*(ampi**2 - p2p3*1D0) + 
     &             p1p2*(ampi**2 - p1p3*1D0 - p2p3*1D0) - 
     &             p1p3*1D0*(ame**2 + p2p3 - p2p4*1D0)) + 
     &          p1p3*2D0*
     &           (-(ampi**2*p1p4*p2p4*1D0) + 
     &             p2p3*(p1p3*(ampi**2 + p2p4) + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0)) + 
     &             p2p3**2*(p1p4 - p1p3*2D0)) + 
     &          kp2*(ampi**2*p1p4*(p2p3 - ampi**2*1D0) + 
     &             ampi**2*p1p3*(p1p4 + p3p4 - p2p3*1D0) + 
     &             p1p3**2*(ampi**2 - p2p3*2D0))) + 
     &       kp3*(-(ame**2*ampi**2*kp4**2*2D0) - 
     &          p1p4*2D0*
     &           (-(ampi**2*p1p4*p2p4*1D0) + 
     &             p2p3*(p1p3*(ampi**2 + p2p4) + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0)) + 
     &             p2p3**2*(p1p4 - p1p3*2D0)) + 
     &          kp4*(ampi**2*p1p4*
     &              (ame**2 + p2p3 - p1p2*1D0 - p2p4*1D0) - 
     &             (ame**2 + p1p2)*p3p4*(ampi**2 - p2p3*1D0)*
     &              2D0 + 
     &             p1p3*(ame**2 + p2p4 - p1p2*1D0 - p2p3*1D0)*
     &              (ampi**2 - p2p3*2D0) - 
     &             ame**2*2D0*
     &              (ampi**2*p2p4 + p2p3*(ampi**2 - p2p3*2D0))
     &             ) + kp2*
     &           (p1p3*(ampi**2*p2p3 - ampi**4*1D0 + 
     &                p1p4*(ampi**2 - p2p3*2D0)) + 
     &             p1p4*(ampi**2*(p1p4 + p2p3) + 
     &                p3p4*(ampi**2 - p2p3*2D0)))) + 
     &       kp1*(-(ampi**2*kp4**2*p2p3*1D0) + 
     &          kp3*(ampi**2*p2p4*(ampi**2 - kp3*1D0) - 
     &             p2p3**2*1D0*(ampi**2 - p3p4*2D0) - 
     &             p2p3*1D0*
     &              (ampi**2*p3p4 + p2p4*(ampi**2 - kp3*2D0))
     &              + p1p4*
     &              (ampi**2*p2p4 + p2p3*(ampi**2 - p2p3*2D0))
     &             ) + kp2*
     &           (ampi**2*kp4*(ampi**2 + p3p4 - p2p3*2D0) + 
     &             kp3*(ampi**4 + p3p4*(ampi**2 - p2p3*2D0)))
     &           + kp4*(ampi**2*p2p3*
     &              (ampi**2 + p1p3 + p2p4 - kp3*3D0) + 
     &             ampi**2*p2p4*(p1p3 - p3p4*1D0 - kp3*3D0) - 
     &             p2p3**2*1D0*(ampi**2 + p1p3*2D0 - kp3*6D0))
     &          )))/(kp1*kp3*kp4*(ame**2 + p1p2)) + 
     &  (C0appx34m2pm1*ep2*2.5D-1*
     &     (kp1**2*(kp2*(ampi**2*kp4 + kp3*p3p4)*2D0 + 
     &          kp4*(-(ampi**2*p2p4*1D0) + 
     &             p2p3*(ampi**2 + p2p3*2D0)) - 
     &          kp3*1D0*(p2p3*(ampi**2 - p2p4*2D0) + 
     &             p2p4*(ampi**2 - p3p4*2D0))) + 
     &       kp2*(ampi**2*kp3**2*(ame**2 + p1p2) + 
     &          kp4*(-(ampi**2*kp2*p1p4*1D0) + 
     &             p1p3**2*p2p4*2D0 + 
     &             kp4*(ame**2*ampi**2 - 
     &                p1p3*(ame**2 + p2p3 - p2p4*1D0)*2D0 + 
     &                p1p2*(ampi**2 - p1p3*2D0)) + 
     &             p1p3*(kp2*
     &                 (-(ampi**2*1D0) + p1p4*2D0 + p3p4*2D0)
     &                 + 2D0*
     &                 (ame**2*(ampi**2 - p3p4*1D0) + 
     &                   p1p2*(ampi**2 - p3p4*1D0) + 
     &                   p1p4*(p2p3 - p2p4*2D0)))) - 
     &          kp3*1D0*(ampi**2*kp2*p1p3 + 
     &             kp4*2D0*
     &              ((ame**2 + p1p2)*p3p4 + 
     &                p1p4*
     &                 (p1p2 + p2p4 - ame**2*1D0 - p2p3*1D0)
     &                 + ame**2*p2p4*2D0) - 
     &             p1p4*1D0*
     &              (p3p4*x12 + kp2*(ampi**2 + p1p4*2D0) - 
     &                2D0*
     &                 (ame**2*ampi**2 + ampi**2*p1p2 + 
     &                   p1p3*p2p4 + p1p4*(p2p3 - p2p4*2D0)))
     &              + ame**2*kp4**2*4D0)) + 
     &       kp1*(kp2**2*(ampi**2*kp3 + kp4*p3p4)*2D0 + 
     &          kp3*2D0*(-(kp4*1D0*
     &                ((ame**2 + p1p2)*p3p4 - 
     &                  p2p3*1D0*(ame**2 + p1p4 - p1p2*1D0) + 
     &                  p1p3*(p2p3 + ame**2*2D0))) + 
     &             p2p4*(p1p4*p2p3 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0) + 
     &                p1p3*(p2p4 - p2p3*2D0))) + 
     &          kp4*(ampi**2*kp4*(ame**2 + p1p2) + 
     &             p2p3*(p3p4*x12 - 
     &                2D0*
     &                 (ame**2*ampi**2 + ampi**2*p1p2 + 
     &                   p1p4*p2p3 + p1p3*(p2p4 - p2p3*2D0))))
     &            + kp3**2*
     &           (ampi**2*(ame**2 + p1p2) - 
     &             p2p4*(ame**2 + p1p2 + p1p4 - p1p3*1D0)*
     &              2D0 - ame**2*kp4*4D0) + 
     &          kp2*(-(kp3**2*p1p4*2D0) + 
     &             kp4*(ampi**2*p2p3 + ampi**2*p2p4 - 
     &                ampi**2*p1p4*1D0 - kp4*p2p3*2D0 - 
     &                p2p4*p3p4*2D0 + 
     &                p1p3*(ampi**2 + p2p3*2D0 + p2p4*2D0)) + 
     &             kp3*(ampi**2*(p1p4 - p2p3*1D0) + 
     &                p2p4*(ampi**2 + p1p4*2D0 - kp4*6D0) + 
     &                p1p3*
     &                 (ampi**2 + p2p4*2D0 - p3p4*2D0 - 
     &                   kp4*6D0))))))/
     &   (kp1*kp2*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D1aex14px34px25m1epm2*ep2*
     &     (ame**2*ampi**2*kp3**2*p1p3 + 
     &       ame**2*ampi**2*kp4**2*p1p3 + 
     &       ame**2*ampi**2*kp3**2*p1p4 + 
     &       ampi**2*kp3**2*p1p2*p1p4 + 
     &       ampi**2*kp4**2*p1p2*p1p4 + 
     &       ame**2*kp3*kp4*p1p4**2 + 
     &       ame**2*ampi**2*kp3**2*p2p4 + 
     &       ame**2*kp3*kp4*p1p3*p2p4 + 
     &       ame**2*kp3*kp4*p1p4*p2p4 + 
     &       ame**2*kp3**2*p2p4**2 - 
     &       ame**2*ampi**2*kp3*kp4*p1p3*1D0 - 
     &       ame**2*kp4**2*p1p3*p1p4*1D0 - 
     &       ame**2*kp3**2*p1p4**2*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*p2p3*1D0 - 
     &       ame**2*kp4**2*p1p3*p2p3*1D0 - 
     &       ame**2*kp4**2*p1p4*p2p3*1D0 - 
     &       ame**2*kp4**2*p2p3**2*1D0 - 
     &       ame**4*ampi**2*kp3*kp4*2D0 + 
     &       ame**2*ampi**2*kp3*kp4*p1p2*2D0 - 
     &       ame**2*kp3*kp4**2*p1p3*2D0 - 
     &       ame**2*kp3**2*kp4*p1p4*2D0 + 
     &       ame**2*kp3*kp4*p1p4*p2p3*2D0 - 
     &       ampi**2*kp4*p1p2*p1p4*p2p3*2D0 - 
     &       kp3*kp4*p1p2*p1p4*p2p3*2D0 + 
     &       ame**2*kp4*p1p3*p1p4*p2p3*2D0 - 
     &       kp3*kp4*p1p3*p1p4*p2p3*2D0 - 
     &       ame**2*kp4*p1p4**2*p2p3*2D0 + 
     &       kp3*kp4*p1p4**2*p2p3*2D0 + 
     &       ame**2*ampi**2*kp4*p2p3**2*2D0 - 
     &       ame**2*kp4*p1p4*p2p3**2*2D0 - 
     &       kp4*p1p4**2*p2p3**2*2D0 - 
     &       ame**2*kp3**2*p1p4*p2p4*2D0 + 
     &       ampi**2*kp3*p1p2*p1p4*p2p4*2D0 - 
     &       kp3**2*p1p2*p1p4*p2p4*2D0 - 
     &       ame**2*kp3*p1p3*p1p4*p2p4*2D0 + 
     &       kp3**2*p1p3*p1p4*p2p4*2D0 + 
     &       ame**2*kp3*p1p4**2*p2p4*2D0 - 
     &       kp3**2*p1p4**2*p2p4*2D0 - 
     &       ame**2*ampi**2*kp3*p2p3*p2p4*2D0 + 
     &       ame**2*kp3*p1p4*p2p3*p2p4*2D0 + 
     &       ame**2*kp4*p1p4*p2p3*p2p4*2D0 - 
     &       kp4*p1p3*p1p4*p2p3*p2p4*2D0 + 
     &       kp3*p1p4**2*p2p3*p2p4*2D0 - 
     &       ame**2*kp3*p1p4*p2p4**2*2D0 + 
     &       kp3*p1p3*p1p4*p2p4**2*2D0 + 
     &       p3p4*(-(ame**2*kp3**2*(p1p4 + p2p4)*1D0) + 
     &          kp4*p2p3*(p1p2*p1p4 - ame**2*p2p3*1D0)*2D0 + 
     &          kp3*(p2p4*(ame**2*p2p3 - p1p2*p1p4*1D0)*2D0 + 
     &             kp4*(-(p1p2*(ame**2 + p1p4)*2D0) + 
     &                ame**2*(p2p3 - p1p3*1D0 + ame**2*2D0))))
     &         + kp1*(kp4*
     &           (ame**2*(ampi**2 - p2p3*1D0)*
     &              (ampi**2 - p3p4*1D0) + 
     &             ame**2*kp3*2D0*
     &              (p3p4 - ampi**2*1D0 + p1p4*2D0) + 
     &             p1p4*(ame**2*p3p4 - 
     &                ampi**2*(ame**2 + p2p4)*1D0 + 
     &                p2p3**2*2D0 + 
     &                p2p3*(ampi**2 + ame**2*2D0))) + 
     &          kp3*(-(ame**2*(p2p4 + p3p4)*1D0*
     &                (ampi**2 - p3p4*1D0)) + 
     &             p1p4*(ame**2*(p3p4 - ampi**2*1D0) - 
     &                p2p3*1D0*(ampi**2 - p2p4*2D0) + 
     &                p2p4*
     &                 (-(ampi**2*1D0) + ame**2*2D0 + 
     &                   p3p4*2D0)))) - 
     &       ame**2*kp3*kp4*p1p3*p1p4*3D0 + 
     &       ame**4*kp3*kp4*p1p4*4D0 - 
     &       ame**2*kp3*kp4*p1p2*p1p4*4D0 + 
     &       kp4*p1p3*p1p4*p2p3**2*4D0 - 
     &       kp3*p1p3*p1p4*p2p3*p2p4*4D0 + 
     &       kp2*(kp3**2*(ame**2*ampi**2 - p1p4**2*2D0) + 
     &          kp4*(-(ame**2*ampi**2*kp4*1D0) - 
     &             ampi**2*p1p4**2*1D0 + 
     &             ame**2*(ampi**2 + p2p3)*
     &              (ampi**2 - p3p4*1D0) + 
     &             p1p4*(ame**2*p3p4 + ampi**2*kp1*2D0 - 
     &                1D0*(ame**2 - p1p3*1D0)*
     &                 (ampi**2 + p2p3*2D0))) + 
     &          kp3*(ampi**2*p1p4**2 + 
     &             ame**2*(ampi**2 - p3p4*1D0)*
     &              (p2p4 - p3p4*1D0) + 
     &             p1p4*(p3p4*(ame**2 + kp1*2D0 - p1p3*2D0) - 
     &                1D0*(ame**2 - p1p3*1D0)*
     &                 (ampi**2 + p2p4*2D0)) + 
     &             kp4*(ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &                p1p4*(ame**2*4D0 - p1p3*6D0))))))/
     &   (kp2*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D1aex13px34px25m1epm2*ep2*
     &     (ame**2*ampi**2*kp4**2*p1p3 + 
     &       ampi**2*kp3**2*p1p2*p1p3 + 
     &       ampi**2*kp4**2*p1p2*p1p3 + 
     &       ame**2*kp3*kp4*p1p3**2 + 
     &       ame**2*ampi**2*kp3**2*p1p4 + 
     &       ame**2*ampi**2*kp4**2*p1p4 + 
     &       ame**2*ampi**2*kp4**2*p2p3 + 
     &       ame**2*kp3*kp4*p1p3*p2p3 + 
     &       ame**2*kp3*kp4*p1p4*p2p3 + 
     &       ame**2*kp4**2*p2p3**2 - 
     &       ame**2*kp4**2*p1p3**2*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*p1p4*1D0 - 
     &       ame**2*kp3**2*p1p3*p1p4*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*p2p4*1D0 - 
     &       ame**2*kp3**2*p1p3*p2p4*1D0 - 
     &       ame**2*kp3**2*p1p4*p2p4*1D0 - 
     &       ame**2*kp3**2*p2p4**2*1D0 - 
     &       ame**4*ampi**2*kp3*kp4*2D0 + 
     &       ame**2*ampi**2*kp3*kp4*p1p2*2D0 - 
     &       ame**2*kp3*kp4**2*p1p3*2D0 - 
     &       ame**2*kp3**2*kp4*p1p4*2D0 - 
     &       ame**2*kp4**2*p1p3*p2p3*2D0 + 
     &       ampi**2*kp4*p1p2*p1p3*p2p3*2D0 - 
     &       kp4**2*p1p2*p1p3*p2p3*2D0 + 
     &       ame**2*kp4*p1p3**2*p2p3*2D0 - 
     &       kp4**2*p1p3**2*p2p3*2D0 - 
     &       ame**2*kp4*p1p3*p1p4*p2p3*2D0 + 
     &       kp4**2*p1p3*p1p4*p2p3*2D0 - 
     &       ame**2*kp4*p1p3*p2p3**2*2D0 + 
     &       kp4*p1p3*p1p4*p2p3**2*2D0 + 
     &       ame**2*kp3*kp4*p1p3*p2p4*2D0 - 
     &       ampi**2*kp3*p1p2*p1p3*p2p4*2D0 - 
     &       kp3*kp4*p1p2*p1p3*p2p4*2D0 - 
     &       ame**2*kp3*p1p3**2*p2p4*2D0 + 
     &       kp3*kp4*p1p3**2*p2p4*2D0 + 
     &       ame**2*kp3*p1p3*p1p4*p2p4*2D0 - 
     &       kp3*kp4*p1p3*p1p4*p2p4*2D0 - 
     &       ame**2*ampi**2*kp4*p2p3*p2p4*2D0 + 
     &       ame**2*kp3*p1p3*p2p3*p2p4*2D0 + 
     &       ame**2*kp4*p1p3*p2p3*p2p4*2D0 + 
     &       kp4*p1p3**2*p2p3*p2p4*2D0 - 
     &       kp3*p1p3*p1p4*p2p3*p2p4*2D0 + 
     &       ame**2*ampi**2*kp3*p2p4**2*2D0 - 
     &       ame**2*kp3*p1p3*p2p4**2*2D0 - 
     &       kp3*p1p3**2*p2p4**2*2D0 + 
     &       p3p4*(-(kp4*1D0*
     &             (ame**2*kp4*(p1p3 + p2p3) + 
     &               p2p3*(p1p2*p1p3 - ame**2*p2p4*1D0)*2D0))
     &           + kp3*(p2p4*(p1p2*p1p3 - ame**2*p2p4*1D0)*
     &              2D0 + 
     &             kp4*(-(p1p2*(ame**2 + p1p3)*2D0) + 
     &                ame**2*(p2p4 - p1p4*1D0 + ame**2*2D0))))
     &         + kp1*(kp3*
     &           (ame**2*(ampi**2 - p2p4*1D0)*
     &              (ampi**2 - p3p4*1D0) + 
     &             ame**2*kp4*2D0*
     &              (p3p4 - ampi**2*1D0 + p1p3*2D0) + 
     &             p1p3*(ampi**2*p2p4 + p2p4**2*2D0 + 
     &                ampi**2*(-(p2p3*1D0) + kp2*2D0) + 
     &                ame**2*(p3p4 - ampi**2*1D0 + p2p4*2D0)))
     &            + kp4*(-(ame**2*(p2p3 + p3p4)*1D0*
     &                (ampi**2 - p3p4*1D0)) + 
     &             p1p3*(-(ampi**2*(ame**2 + p2p4)*1D0) + 
     &                p3p4*(ame**2 + kp2*2D0) + 
     &                p2p3*
     &                 (-(ampi**2*1D0) + ame**2*2D0 + 
     &                   p2p4*2D0 + p3p4*2D0)))) - 
     &       ame**2*kp3*kp4*p1p3*p1p4*3D0 + 
     &       ame**4*kp3*kp4*p1p3*4D0 - 
     &       ame**2*kp3*kp4*p1p2*p1p3*4D0 - 
     &       kp4*p1p3*p1p4*p2p3*p2p4*4D0 + 
     &       kp3*p1p3*p1p4*p2p4**2*4D0 + 
     &       kp2*(-(ame**2*ampi**2*kp3**2*1D0) + 
     &          kp4*(ampi**2*p1p3**2 + 
     &             ame**2*(ampi**2 - p3p4*1D0)*
     &              (p2p3 - p3p4*1D0) + 
     &             kp4*(ame**2*ampi**2 - p1p3**2*2D0) + 
     &             p1p3*(p3p4*(ame**2 - p1p4*2D0) - 
     &                1D0*(ame**2 - p1p4*1D0)*
     &                 (ampi**2 + p2p3*2D0))) + 
     &          kp3*(-(ampi**2*p1p3**2*1D0) + 
     &             ame**2*(ampi**2 + p2p4)*
     &              (ampi**2 - p3p4*1D0) + 
     &             p1p3*(ame**2*p3p4 - 
     &                1D0*(ame**2 - p1p4*1D0)*
     &                 (ampi**2 + p2p4*2D0)) + 
     &             kp4*(ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &                p1p3*(ame**2*4D0 - p1p4*6D0))))))/
     &   (kp2*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D1aex24px34px15m2epm1*ep2*
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
     &       kp1*(kp3**2*(ame**2*ampi**2 - p2p4**2*2D0) - 
     &          kp4*1D0*(ame**2*ampi**2*kp4 + 
     &             ampi**2*p2p4**2 + ame**2*ampi**2*p3p4 + 
     &             ampi**2*p2p4*x25 - ame**2*ampi**4*1D0 - 
     &             ampi**2*p2p3*p2p4*1D0 - 
     &             ame**2*p2p4*p3p4*1D0 + 
     &             p1p3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &                p2p4*(ame**2 - p2p3*1D0)*2D0)) + 
     &          kp3*(ame**2*p3p4**2 + 
     &             ampi**2*p2p4*(p2p3 + p2p4 - ame**2*1D0) + 
     &             p1p4*(ame**2*(ampi**2 - p3p4*1D0) + 
     &                p2p4*(p2p3 - ame**2*1D0)*2D0) + 
     &             p3p4*(-(ame**2*ampi**2*1D0) + 
     &                p2p4*(ame**2 + kp2*2D0 - p2p3*2D0)) + 
     &             kp4*(ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &                p2p4*(ame**2*4D0 - p2p3*6D0))))))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) + 
     &  (D1aex23px34px15m2epm1*ep2*
     &     (ame**2*ampi**2*kp4**2*p1p3 + 
     &       ame**2*kp4**2*p1p3**2 + 
     &       ame**2*ampi**2*kp4**2*p2p3 + 
     &       ampi**2*kp3**2*p1p2*p2p3 + 
     &       ampi**2*kp4**2*p1p2*p2p3 + 
     &       ame**2*kp3*kp4*p1p3*p2p3 + 
     &       ame**2*kp3*kp4*p2p3**2 + 
     &       ame**2*ampi**2*kp3**2*p2p4 + 
     &       ame**2*ampi**2*kp4**2*p2p4 + 
     &       ame**2*kp3*kp4*p1p3*p2p4 - 
     &       ame**2*ampi**2*kp3*kp4*p1p4*1D0 - 
     &       ame**2*kp3**2*p1p4**2*1D0 - 
     &       ame**2*kp3**2*p1p4*p2p3*1D0 - 
     &       ame**2*kp4**2*p2p3**2*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*p2p4*1D0 - 
     &       ame**2*kp3**2*p1p4*p2p4*1D0 - 
     &       ame**2*kp3**2*p2p3*p2p4*1D0 - 
     &       ame**4*ampi**2*kp3*kp4*2D0 + 
     &       ame**2*ampi**2*kp3*kp4*p1p2*2D0 - 
     &       ame**2*ampi**2*kp4*p1p3*p1p4*2D0 + 
     &       ame**2*ampi**2*kp3*p1p4**2*2D0 - 
     &       ame**2*kp3*kp4**2*p2p3*2D0 - 
     &       ame**2*kp4**2*p1p3*p2p3*2D0 + 
     &       ampi**2*kp4*p1p2*p1p3*p2p3*2D0 - 
     &       kp4**2*p1p2*p1p3*p2p3*2D0 - 
     &       ame**2*kp4*p1p3**2*p2p3*2D0 + 
     &       ame**2*kp3*kp4*p1p4*p2p3*2D0 - 
     &       ampi**2*kp3*p1p2*p1p4*p2p3*2D0 - 
     &       kp3*kp4*p1p2*p1p4*p2p3*2D0 + 
     &       ame**2*kp3*p1p3*p1p4*p2p3*2D0 + 
     &       ame**2*kp4*p1p3*p1p4*p2p3*2D0 - 
     &       ame**2*kp3*p1p4**2*p2p3*2D0 + 
     &       ame**2*kp4*p1p3*p2p3**2*2D0 - 
     &       kp4**2*p1p3*p2p3**2*2D0 - 
     &       ame**2*kp3*p1p4*p2p3**2*2D0 + 
     &       kp3*kp4*p1p4*p2p3**2*2D0 + 
     &       kp4*p1p3*p1p4*p2p3**2*2D0 - 
     &       kp3*p1p4**2*p2p3**2*2D0 - 
     &       ame**2*kp3**2*kp4*p2p4*2D0 - 
     &       ame**2*kp4*p1p3*p2p3*p2p4*2D0 + 
     &       kp4**2*p1p3*p2p3*p2p4*2D0 + 
     &       kp4*p1p3**2*p2p3*p2p4*2D0 + 
     &       ame**2*kp3*p1p4*p2p3*p2p4*2D0 - 
     &       kp3*kp4*p1p4*p2p3*p2p4*2D0 - 
     &       kp3*p1p3*p1p4*p2p3*p2p4*2D0 + 
     &       p3p4*(kp4*(-(ame**2*kp4*(p1p3 + p2p3)*1D0) + 
     &             p1p3*(ame**2*p1p4 - p1p2*p2p3*1D0)*2D0) + 
     &          kp3*(p1p4*(p1p2*p2p3 - ame**2*p1p4*1D0)*2D0 + 
     &             kp4*(-(p1p2*(ame**2 + p2p3)*2D0) + 
     &                ame**2*(p1p4 - p2p4*1D0 + ame**2*2D0))))
     &         + kp2*(kp3*
     &           (ame**2*ampi**4 - 
     &             ampi**2*(ame**2 + p1p3)*p2p3*1D0 + 
     &             ame**2*p3p4*(p2p3 - ampi**2*1D0) + 
     &             p1p4**2*p2p3*2D0 + 
     &             ame**2*kp4*2D0*
     &              (p3p4 - ampi**2*1D0 + p2p3*2D0) + 
     &             p1p4*(ame**2*(p3p4 - ampi**2*1D0) + 
     &                p2p3*(ampi**2 + ame**2*2D0))) + 
     &          kp4*(ame**2*p3p4*(p3p4 - ampi**2*1D0) + 
     &             p2p3*(ame**2*p3p4 - 
     &                ampi**2*(ame**2 + p1p4)*1D0) + 
     &             p1p3*(ame**2*(p3p4 - ampi**2*1D0) + 
     &                p2p3*
     &                 (-(ampi**2*1D0) + ame**2*2D0 + 
     &                   p1p4*2D0 + p3p4*2D0)))) - 
     &       ame**2*kp3*kp4*p2p3*p2p4*3D0 + 
     &       ame**4*kp3*kp4*p2p3*4D0 - 
     &       ame**2*kp3*kp4*p1p2*p2p3*4D0 - 
     &       kp4*p1p3*p1p4*p2p3*p2p4*4D0 + 
     &       kp3*p1p4**2*p2p3*p2p4*4D0 + 
     &       kp1*(-(ame**2*ampi**2*kp3**2*1D0) + 
     &          kp4*(ame**2*p3p4**2 + 
     &             ampi**2*p2p3*(p2p3 + p2p4 - ame**2*1D0) + 
     &             kp4*(ame**2*ampi**2 - p2p3**2*2D0) + 
     &             p1p3*(ame**2*(ampi**2 - p3p4*1D0) + 
     &                p2p3*(p2p4 - ame**2*1D0)*2D0) + 
     &             p3p4*(-(ame**2*ampi**2*1D0) + 
     &                p2p3*(ame**2 + kp2*2D0 - p2p4*2D0))) + 
     &          kp3*(ame**2*ampi**4 + ampi**2*p2p3*p2p4 + 
     &             ame**2*p2p3*p3p4 - ampi**2*p2p3**2*1D0 - 
     &             ame**2*ampi**2*p3p4*1D0 - 
     &             ampi**2*p2p3*x25*1D0 + 
     &             p1p4*(ame**2*(ampi**2 - p3p4*1D0) + 
     &                p2p3*(p2p4 - ame**2*1D0)*2D0) + 
     &             kp4*(ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &                p2p3*(ame**2*4D0 - p2p4*6D0))))))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)) - 
     &  (ame**2*D3aex15ex250x34em1em2*ep2*2D0*
     &     (kp3*kp4*p1p3 + kp3*kp4*p1p4 + kp3*kp4*p2p3 - 
     &       kp4**2*p1p3*1D0 - kp3**2*p1p4*1D0 - 
     &       kp4**2*p2p3*1D0 + kp4*p1p3*p2p3*4D0 - 
     &       kp4*p1p4*p2p3*4D0 - kp4*p2p3**2*4D0 - 
     &       kp3*p2p4**2*4D0 + 
     &       p2p4*((kp3 + kp4)*p2p3*4D0 + 
     &          kp3*(kp4 - kp3*1D0 - p1p3*4D0 + p1p4*4D0)) + 
     &       kp2*(kp4*(p3p4 - ampi**2*1D0 - p2p3*4D0) + 
     &          kp3*(p3p4 - ampi**2*1D0 + kp4*4D0 - p2p4*4D0))
     &         + kp1*(kp4*(p3p4 - ampi**2*1D0 + p2p3*4D0) + 
     &          kp3*(p3p4 - ampi**2*1D0 + kp4*4D0 + p2p4*4D0))
     &         + ame**2*kp3*kp4*8D0 - kp3*kp4*p1p2*8D0))/
     &   (kp3*kp4*(ame**2 + p1p2)) - 
     &  (D1aex15ex250x34em1em2*ep2*2D0*
     &     (ame**2*kp4**2*p1p3 + ame**2*kp3**2*p1p4 + 
     &       ame**2*kp3*kp4*p2p3 - ame**2*kp3*kp4*p1p3*1D0 - 
     &       ame**2*kp3*kp4*p1p4*1D0 - 
     &       ame**2*kp4**2*p2p3*1D0 - kp3*kp4*p1p2*p1p3*2D0 + 
     &       kp4**2*p1p2*p1p3*2D0 + kp3**2*p1p2*p1p4*2D0 - 
     &       kp3*kp4*p1p2*p1p4*2D0 - 
     &       kp1**2*(kp4*p2p3 + kp3*p2p4)*2D0 + 
     &       kp2*(kp3*(ame**2*p3p4 - ame**2*ampi**2*1D0 + 
     &             p1p4*2D0*(p1p4 - p1p3*1D0 + ame**2*2D0)) + 
     &          kp4*(ame**2*p3p4 - ame**2*ampi**2*1D0 + 
     &             p1p3*2D0*(p1p3 - p1p4*1D0 + ame**2*2D0)))
     &        - ame**2*kp4*p1p3**2*4D0 + 
     &       ame**2*kp3*p1p3*p1p4*4D0 + 
     &       ame**2*kp4*p1p3*p1p4*4D0 - 
     &       ame**2*kp3*p1p4**2*4D0 + 
     &       ame**2*kp4*p1p3*p2p3*4D0 - 
     &       ame**2*kp3*p1p4*p2p3*4D0 + 
     &       kp1*(kp4*(ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p4*p2p3*2D0 + 
     &             p1p2*(ampi**2 - p3p4*1D0)*2D0 - 
     &             p1p3*2D0*(p2p3 - kp2*1D0 + ame**2*2D0)) + 
     &          kp3*(p1p3*p2p4*2D0 - 
     &             p1p4*2D0*(p2p4 - kp2*1D0 + ame**2*2D0) + 
     &             ampi**2*(ame**2 + p1p2*2D0) - 
     &             p3p4*1D0*(ame**2 + p1p2*2D0) - 
     &             kp4*(ame**2 + p1p2)*4D0)) + 
     &       ame**2*p2p4*
     &        (-(kp3**2*1D0) - kp4*p1p3*4D0 + 
     &          kp3*(kp4 + p1p4*4D0)) + ame**4*kp3*kp4*8D0 - 
     &       ame**2*kp3*kp4*p1p2*8D0))/
     &   (kp3*kp4*(ame**2 + p1p2)) + 
     &  (D1aeex35px12x14m1em2p*ep3*2D0*
     &     (ame**2*ampi**2*kp3**2*p1p3 + 
     &       ame**2*ampi**2*kp3**2*p1p4 + 
     &       ame**2*ampi**2*kp3*kp4*p1p4 + 
     &       ampi**2*kp4**2*p1p2*p1p4 + 
     &       ame**2*ampi**2*kp4*p1p3*p1p4 + 
     &       ame**2*ampi**4*kp4*p2p3 + 
     &       ame**2*ampi**2*kp3**2*p2p4 + 
     &       ame**2*ampi**2*kp4*p1p4*p2p4 - 
     &       ampi**4*kp4*p1p2*p1p4*1D0 - 
     &       ame**2*ampi**2*kp4*p1p4**2*1D0 - 
     &       ame**2*ampi**2*kp3**2*p2p3*1D0 - 
     &       ame**2*ampi**2*kp4**2*p2p3*1D0 - 
     &       ame**2*ampi**2*kp4*p1p4*p2p3*1D0 - 
     &       ampi**2*kp4*p1p4**2*p2p3*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*p2p4*1D0 - 
     &       ampi**2*kp4*p1p3*p1p4*p2p4*1D0 + 
     &       kp3*p3p4**2*(ame**2*p2p3 - p1p2*p1p4*1D0) + 
     &       ampi**2*kp3**2*p1p2*p1p4*2D0 + 
     &       ampi**2*kp3*kp4*p1p2*p1p4*2D0 - 
     &       kp3**2*kp4*p1p2*p1p4*2D0 - 
     &       ame**2*kp3*kp4*p1p3*p1p4*2D0 - 
     &       ame**2*ampi**2*kp3*kp4*p2p3*2D0 + 
     &       ame**2*kp3*kp4**2*p2p3*2D0 + 
     &       ame**2*kp3*kp4*p1p4*p2p3*2D0 + 
     &       ampi**2*kp4*p1p3*p1p4*p2p3*2D0 + 
     &       kp3**2*p1p4**2*p2p3*2D0 + 
     &       kp3*kp4*p1p4**2*p2p3*2D0 + 
     &       ame**2*kp3**2*kp4*p2p4*2D0 + 
     &       kp3**2*p1p3*p1p4*p2p4*2D0 + 
     &       kp3*kp4*p1p3*p1p4*p2p4*2D0 - 
     &       ame**2*kp3**2*kp4*p1p4*4D0 - 
     &       kp3**2*p1p4**2*p2p4*4D0 + 
     &       kp1*(ampi**2*kp4*p1p4*
     &           (ame**2 - p2p4*1D0 + kp2*2D0 + p2p3*2D0) - 
     &          kp3**2*1D0*(ame**2*ampi**2 + p1p4*p2p4*2D0) + 
     &          kp3*(p3p4*
     &              (ame**2*(p3p4 - ampi**2*1D0) + 
     &                p1p4*(ame**2 + p2p4 + kp2*2D0)) - 
     &             kp4*p1p4*p2p3*4D0)) + 
     &       kp2*(kp3**2*(ame**2*ampi**2 - p1p4**2*2D0) + 
     &          ampi**2*kp4*
     &           (ame**2*ampi**2 - ame**2*kp4*1D0 - 
     &             ame**2*p3p4*1D0 - 
     &             p1p4*1D0*(ame**2 + p1p4 - p1p3*2D0)) + 
     &          kp3*(p1p4*p3p4*(p1p4 - ame**2*1D0) + 
     &             kp4*(ame**2*p3p4 - p1p3*p1p4*4D0))) + 
     &       p3p4*(ampi**2*kp4*
     &           (p1p2*p1p4 - ame**2*p2p3*1D0) - 
     &          kp3**2*(ame**2 + p1p2)*p1p4*2D0 + 
     &          kp3*(ampi**2*p1p2*p1p4 + ame**2*p1p4**2 + 
     &             ame**2*p1p4*p2p3 + p1p4**2*p2p3 - 
     &             ame**2*ampi**2*p2p3*1D0 - 
     &             ame**2*p1p4*p2p4*1D0 + 
     &             p1p3*p1p4*(p2p4 - ame**2*1D0 - p2p3*2D0) - 
     &             kp4*1D0*
     &              (p1p2*p1p4*3D0 + ame**2*(p1p3 - p2p3*4D0))
     &             )) - kp3*kp4*p1p3*p1p4*p2p3*8D0))/
     &   (kp3**2*kp4*(ame**2 + p1p2)) - 
     &  (D1aeex45px12x13m1em2p*ep3*2D0*
     &     (ame**2*ampi**2*kp3*kp4*p1p3 + 
     &       ame**2*ampi**2*kp4**2*p1p3 + 
     &       ampi**2*kp3**2*p1p2*p1p3 + 
     &       ame**2*ampi**2*kp4**2*p1p4 + 
     &       ame**2*ampi**2*kp3*p1p3*p1p4 + 
     &       ame**2*ampi**2*kp4**2*p2p3 + 
     &       ame**2*ampi**2*kp3*p1p3*p2p3 + 
     &       ame**2*ampi**4*kp3*p2p4 - 
     &       ampi**4*kp3*p1p2*p1p3*1D0 - 
     &       ame**2*ampi**2*kp3*p1p3**2*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*p2p3*1D0 - 
     &       ampi**2*kp3*p1p3*p1p4*p2p3*1D0 - 
     &       ame**2*ampi**2*kp3**2*p2p4*1D0 - 
     &       ame**2*ampi**2*kp4**2*p2p4*1D0 - 
     &       ame**2*ampi**2*kp3*p1p3*p2p4*1D0 - 
     &       ampi**2*kp3*p1p3**2*p2p4*1D0 + 
     &       kp4*p3p4**2*(ame**2*p2p4 - p1p2*p1p3*1D0) + 
     &       ampi**2*kp3*kp4*p1p2*p1p3*2D0 + 
     &       ampi**2*kp4**2*p1p2*p1p3*2D0 - 
     &       kp3*kp4**2*p1p2*p1p3*2D0 - 
     &       ame**2*kp3*kp4*p1p3*p1p4*2D0 + 
     &       ame**2*kp3*kp4**2*p2p3*2D0 + 
     &       kp3*kp4*p1p3*p1p4*p2p3*2D0 + 
     &       kp4**2*p1p3*p1p4*p2p3*2D0 - 
     &       ame**2*ampi**2*kp3*kp4*p2p4*2D0 + 
     &       ame**2*kp3**2*kp4*p2p4*2D0 + 
     &       ame**2*kp3*kp4*p1p3*p2p4*2D0 + 
     &       kp3*kp4*p1p3**2*p2p4*2D0 + 
     &       kp4**2*p1p3**2*p2p4*2D0 + 
     &       ampi**2*kp3*p1p3*p1p4*p2p4*2D0 + 
     &       kp2*(-(ame**2*ampi**2*kp3**2*1D0) + 
     &          kp4*(p1p3*p3p4*(p1p3 - x15*1D0) + 
     &             kp4*(ame**2*ampi**2 - p1p3**2*2D0)) + 
     &          kp3*(-(ampi**2*p1p3**2*1D0) + 
     &             ame**2*
     &              (ampi**4 + p3p4*(kp4 - ampi**2*1D0)) + 
     &             p1p3*(-(ampi**2*x15*1D0) + 
     &                p1p4*2D0*(ampi**2 - kp4*2D0)))) - 
     &       ame**2*kp3*kp4**2*p1p3*4D0 - 
     &       kp4**2*p1p3**2*p2p3*4D0 + 
     &       kp1*(-(kp4**2*1D0*
     &             (ame**2*ampi**2 + p1p3*p2p3*2D0)) + 
     &          ampi**2*kp3*p1p3*
     &           (ame**2 - p2p3*1D0 + p2p4*2D0) + 
     &          kp4*(p3p4*
     &              (p1p3*(ame**2 + p2p3) + 
     &                ame**2*(p3p4 - ampi**2*1D0)) - 
     &             kp3*p1p3*p2p4*4D0)) + 
     &       p3p4*(ampi**2*kp3*
     &           (p1p2*p1p3 - ame**2*p2p4*1D0) - 
     &          kp4**2*(ame**2 + p1p2)*p1p3*2D0 + 
     &          kp4*(p1p3**2*(ame**2 + p2p4) + 
     &             p1p3*(ame**2*(p2p4 - p2p3*1D0) + 
     &                p1p4*(p2p3 - ame**2*1D0 - p2p4*2D0)) + 
     &             p1p2*p1p3*(ampi**2 - kp3*3D0) - 
     &             ame**2*1D0*
     &              (ampi**2*p2p4 + kp3*(p1p4 - p2p4*4D0))))
     &        - kp3*kp4*p1p3*p1p4*p2p4*8D0))/
     &   (kp3*kp4**2*(ame**2 + p1p2)) - 
     &  (D1aeepx35x12x24m1em2p*ep3*2D0*
     &     (ame**2*ampi**4*kp4*p1p3 + 
     &       ame**2*ampi**2*kp3**2*p1p4 + 
     &       ame**2*ampi**2*kp3**2*p2p3 + 
     &       ame**2*ampi**2*kp3**2*p2p4 + 
     &       ame**2*ampi**2*kp3*kp4*p2p4 + 
     &       ampi**2*kp4**2*p1p2*p2p4 + 
     &       ame**2*ampi**2*kp4*p1p4*p2p4 + 
     &       ame**2*ampi**2*kp4*p2p3*p2p4 + 
     &       ampi**2*kp3*p1p2*p2p4*p3p4 + 
     &       ampi**2*kp4*p1p2*p2p4*p3p4 + 
     &       ame**2*kp3*p1p3*p2p4*p3p4 + 
     &       kp3*p1p4*p2p3*p2p4*p3p4 + 
     &       ame**2*kp3*p2p4**2*p3p4 + 
     &       kp3*p1p3*p2p4**2*p3p4 + 
     &       ame**2*kp3*p1p3*p3p4**2 - 
     &       ame**2*ampi**2*kp3**2*p1p3*1D0 - 
     &       ame**2*ampi**2*kp4**2*p1p3*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*p1p4*1D0 - 
     &       ampi**4*kp4*p1p2*p2p4*1D0 - 
     &       ampi**2*kp4*p1p4*p2p3*p2p4*1D0 - 
     &       ame**2*ampi**2*kp4*p2p4**2*1D0 - 
     &       ampi**2*kp4*p1p3*p2p4**2*1D0 - 
     &       ame**2*ampi**2*kp3*p1p3*p3p4*1D0 - 
     &       ame**2*ampi**2*kp4*p1p3*p3p4*1D0 - 
     &       ame**2*kp3*kp4*p2p3*p3p4*1D0 - 
     &       ame**2*kp3*p1p4*p2p4*p3p4*1D0 - 
     &       ame**2*kp3*p2p3*p2p4*p3p4*1D0 - 
     &       kp3*p1p2*p2p4*p3p4**2*1D0 - 
     &       ampi**2*kp4*p1p3*p2p4*x25*1D0 - 
     &       ame**2*ampi**2*kp3*kp4*p1p3*2D0 + 
     &       ame**2*kp3*kp4**2*p1p3*2D0 + 
     &       ame**2*kp3**2*kp4*p1p4*2D0 + 
     &       ampi**2*kp3**2*p1p2*p2p4*2D0 + 
     &       ampi**2*kp3*kp4*p1p2*p2p4*2D0 - 
     &       kp3**2*kp4*p1p2*p2p4*2D0 + 
     &       ame**2*kp3*kp4*p1p3*p2p4*2D0 - 
     &       ame**2*kp3*kp4*p2p3*p2p4*2D0 + 
     &       ampi**2*kp4*p1p3*p2p3*p2p4*2D0 + 
     &       kp3**2*p1p4*p2p3*p2p4*2D0 + 
     &       kp3*kp4*p1p4*p2p3*p2p4*2D0 + 
     &       kp3**2*p1p3*p2p4**2*2D0 + 
     &       kp3*kp4*p1p3*p2p4**2*2D0 - 
     &       ame**2*kp3**2*p2p4*p3p4*2D0 - 
     &       kp3**2*p1p2*p2p4*p3p4*2D0 - 
     &       kp3*p1p3*p2p3*p2p4*p3p4*2D0 - 
     &       kp3*kp4*p1p2*p2p4*p3p4*3D0 - 
     &       ame**2*kp3**2*kp4*p2p4*4D0 - 
     &       kp3**2*p1p4*p2p4**2*4D0 + 
     &       ame**2*kp3*kp4*p1p3*p3p4*4D0 + 
     &       kp2*(ampi**2*kp4*p2p4*(ame**2 - p1p4*1D0) - 
     &          kp3**2*1D0*(ame**2*ampi**2 + p1p4*p2p4*2D0) + 
     &          kp3*(p3p4*
     &              ((ame**2 + p1p4)*p2p4 + 
     &                ame**2*(p3p4 - ampi**2*1D0)) - 
     &             kp4*p1p3*p2p4*4D0)) + 
     &       kp1*(kp3**2*(ame**2*ampi**2 - p2p4**2*2D0) + 
     &          ampi**2*kp4*
     &           (ame**2*ampi**2 - ame**2*kp4*1D0 - 
     &             ame**2*p3p4*1D0 - 
     &             p2p4*1D0*
     &              (ame**2 + p2p4 - kp2*2D0 - p2p3*2D0)) + 
     &          kp3*(p2p4*p3p4*
     &              (p2p4 - ame**2*1D0 + kp2*2D0) + 
     &             kp4*(ame**2*p3p4 - p2p3*p2p4*4D0))) - 
     &       kp3*kp4*p1p3*p2p3*p2p4*8D0))/
     &   (kp3**2*kp4*(ame**2 + p1p2)) + 
     &  (D1aeepx45x12x23m1em2p*ep3*2D0*
     &     (ame**2*ampi**2*kp4**2*p1p3 + 
     &       ame**2*ampi**4*kp3*p1p4 + 
     &       ame**2*ampi**2*kp3*kp4*p2p3 + 
     &       ame**2*ampi**2*kp4**2*p2p3 + 
     &       ampi**2*kp3**2*p1p2*p2p3 + 
     &       ame**2*ampi**2*kp3*p1p3*p2p3 + 
     &       ame**2*ampi**2*kp4**2*p2p4 + 
     &       ame**2*ampi**2*kp3*p2p3*p2p4 + 
     &       ampi**2*kp3*p1p2*p2p3*p3p4 + 
     &       ampi**2*kp4*p1p2*p2p3*p3p4 + 
     &       ame**2*kp4*p1p4*p2p3*p3p4 + 
     &       ame**2*kp4*p2p3**2*p3p4 + 
     &       kp4*p1p4*p2p3**2*p3p4 + 
     &       kp4*p1p3*p2p3*p2p4*p3p4 + 
     &       ame**2*kp4*p1p4*p3p4**2 - 
     &       ame**2*ampi**2*kp3*kp4*p1p3*1D0 - 
     &       ame**2*ampi**2*kp3**2*p1p4*1D0 - 
     &       ame**2*ampi**2*kp4**2*p1p4*1D0 - 
     &       ampi**4*kp3*p1p2*p2p3*1D0 - 
     &       ame**2*ampi**2*kp3*p2p3**2*1D0 - 
     &       ampi**2*kp3*p1p4*p2p3**2*1D0 - 
     &       ampi**2*kp3*p1p3*p2p3*p2p4*1D0 - 
     &       ame**2*ampi**2*kp3*p1p4*p3p4*1D0 - 
     &       ame**2*ampi**2*kp4*p1p4*p3p4*1D0 - 
     &       ame**2*kp4*p1p3*p2p3*p3p4*1D0 - 
     &       ame**2*kp3*kp4*p2p4*p3p4*1D0 - 
     &       ame**2*kp4*p2p3*p2p4*p3p4*1D0 - 
     &       kp4*p1p2*p2p3*p3p4**2*1D0 - 
     &       ampi**2*kp3*p1p4*p2p3*x25*1D0 + 
     &       ame**2*kp3*kp4**2*p1p3*2D0 - 
     &       ame**2*ampi**2*kp3*kp4*p1p4*2D0 + 
     &       ame**2*kp3**2*kp4*p1p4*2D0 + 
     &       ampi**2*kp3*kp4*p1p2*p2p3*2D0 + 
     &       ampi**2*kp4**2*p1p2*p2p3*2D0 - 
     &       kp3*kp4**2*p1p2*p2p3*2D0 + 
     &       ame**2*kp3*kp4*p1p4*p2p3*2D0 + 
     &       kp3*kp4*p1p4*p2p3**2*2D0 + 
     &       kp4**2*p1p4*p2p3**2*2D0 - 
     &       ame**2*kp3*kp4*p2p3*p2p4*2D0 + 
     &       kp3*kp4*p1p3*p2p3*p2p4*2D0 + 
     &       kp4**2*p1p3*p2p3*p2p4*2D0 + 
     &       ampi**2*kp3*p1p4*p2p3*p2p4*2D0 - 
     &       ame**2*kp4**2*p2p3*p3p4*2D0 - 
     &       kp4**2*p1p2*p2p3*p3p4*2D0 - 
     &       kp4*p1p4*p2p3*p2p4*p3p4*2D0 + 
     &       kp1*(-(ame**2*ampi**2*kp3**2*1D0) + 
     &          kp4*(p2p3*p3p4*
     &              (p2p3 - ame**2*1D0 + kp2*2D0) + 
     &             kp4*(ame**2*ampi**2 - p2p3**2*2D0)) + 
     &          kp3*(-(ampi**2*p2p3**2*1D0) + 
     &             ame**2*
     &              (ampi**4 + p3p4*(kp4 - ampi**2*1D0)) + 
     &             p2p3*(-(ampi**2*x25*1D0) + 
     &                p2p4*2D0*(ampi**2 - kp4*2D0)))) - 
     &       kp3*kp4*p1p2*p2p3*p3p4*3D0 - 
     &       ame**2*kp3*kp4**2*p2p3*4D0 - 
     &       kp4**2*p1p3*p2p3**2*4D0 + 
     &       ame**2*kp3*kp4*p1p4*p3p4*4D0 + 
     &       kp2*(ampi**2*kp3*p2p3*(ame**2 - p1p3*1D0) - 
     &          kp4**2*1D0*(ame**2*ampi**2 + p1p3*p2p3*2D0) + 
     &          kp4*(p3p4*
     &              ((ame**2 + p1p3)*p2p3 + 
     &                ame**2*(p3p4 - ampi**2*1D0)) - 
     &             kp3*p1p4*p2p3*4D0)) - 
     &       kp3*kp4*p1p4*p2p3*p2p4*8D0))/
     &   (kp3*kp4**2*(ame**2 + p1p2)) + 
     &  (D2aex12x35x14epem1m2p*ep3*2D0*
     &     (kp2**2*p1p4*(ampi**2*kp4 + kp3*p3p4) + 
     &       ampi**2*kp4*
     &        (-(p1p4**2*p2p3*1D0) + 
     &          p2p3*(ampi**2*(ame**2 + p1p2) + 
     &             (kp1 + p1p3)*p2p4 - 
     &             (ame**2 + p1p2)*(kp4 + p3p4)*1D0) + 
     &          p1p4*((ame**2 + p1p2)*(kp4 + p3p4) - 
     &             ame**2*ampi**2*1D0 - ampi**2*p1p2*1D0 + 
     &             p2p3*(p2p3 + (kp1 + p1p3)*2D0) - 
     &             p2p4*1D0*(kp1 + p1p3 + p2p3*2D0))) + 
     &       kp3*(-((ame**2 + p1p2)*p3p4**2*1D0*
     &             (p1p4 - p2p3*1D0)) + 
     &          kp4*2D0*(p1p4**2*p2p3 + 
     &             p2p3*(-((ame**2 + p1p2)*1D0*
     &                   (ampi**2 - kp4*1D0)) - 
     &                p2p4*1D0*(kp1 + p1p3*2D0)) + 
     &             p1p4*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &                p2p3*2D0*(kp1 + p2p3 + p1p3*2D0) + 
     &                p2p4*(p1p3 + p2p3*2D0))) + 
     &          p3p4*(p1p4**2*p2p3 + 
     &             p1p4*(ame**2*ampi**2 + ampi**2*p1p2 - 
     &                p2p3**2*1D0 - p1p3*p2p3*2D0 + 
     &                p2p4*(kp1 + p1p3 + p2p3*2D0) - 
     &                kp4*(ame**2 + p1p2)*3D0) + 
     &             p2p3*(p2p4*(kp1 - p1p3*1D0) - 
     &                (ame**2 + p1p2)*1D0*(ampi**2 - kp4*3D0))
     &             )) + kp2*
     &        (-(ampi**2*kp4**2*(ame**2 + p1p2)*1D0) + 
     &          kp3*(-(kp3*p1p4*(p1p4 + p2p4)*2D0) + 
     &             p3p4*(ame**2*(p3p4 - ampi**2*1D0) + 
     &                p1p2*(p3p4 - ampi**2*1D0) + 
     &                p1p4*(p1p4 + kp1*2D0) + 
     &                p2p4*(kp1 - p1p3*1D0 + p1p4*2D0))) + 
     &          kp4*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &             ampi**2*kp1*p2p4 + ampi**2*p1p3*p2p4 - 
     &             ampi**2*p1p4**2*1D0 - 
     &             (ame**2 + p1p2)*p3p4*1D0*
     &              (ampi**2 - kp3*1D0) + 
     &             ampi**2*kp1*p1p4*2D0 + 
     &             ampi**2*p1p3*p1p4*2D0 + 
     &             ampi**2*p1p4*p2p3*2D0 - 
     &             kp3*p1p3*p2p4*2D0 - 
     &             ampi**2*p1p4*p2p4*2D0 - 
     &             kp3*p1p3*p1p4*4D0 - kp3*p1p4*p2p3*4D0)) + 
     &       kp3**2*(p2p4*
     &           (ampi**2*(ame**2 + p1p2) - 
     &             (ame**2 + p1p2)*p3p4*1D0 + 
     &             kp4*(p1p2 - ame**2*1D0) + 
     &             p2p4*(p1p3 - kp1*1D0))*2D0 + 
     &          p1p4**2*2D0*(p2p3 - p2p4*2D0) + 
     &          p1p4*(ampi**2*(ame**2 + p1p2)*2D0 - 
     &             (ame**2 + p1p2)*p3p4*2D0 + 
     &             p2p4*(p1p3 + p2p3 - kp1*1D0)*2D0 - 
     &             p2p4**2*4D0 + 
     &             kp4*(x12 - p1p2*4D0 - ame**2*8D0)))))/
     &   (kp3**2*kp4*(ame**2 + p1p2)) - 
     &  (D2aex12x45x13epem1m2p*ep3*2D0*
     &     (kp2**2*p1p3*(ampi**2*kp3 + kp4*p3p4) + 
     &       kp3**2*(ame**2 + p1p2)*
     &        (ampi**2*p1p3 - p2p4*1D0*(ampi**2 - kp4*2D0)) + 
     &       kp4*(kp4*(p1p3 + p2p3)*2D0*
     &           (p1p3*p2p4 + ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0) + 
     &             p2p3*(p1p4 - kp1*1D0 - p1p3*2D0)) + 
     &          p3p4*(p1p3**2*p2p4 + 
     &             p2p4*(ame**2*(p3p4 - ampi**2*1D0) + 
     &                p1p2*(p3p4 - ampi**2*1D0) + 
     &                p2p3*(kp1 - p1p4*1D0)) + 
     &             p1p3*((kp1 + p1p4)*p2p3 - p2p4**2*1D0 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0) + 
     &                p2p4*(p2p3 - p1p4*1D0)*2D0))) + 
     &       kp2*(-(ampi**2*kp3**2*(ame**2 + p1p2)*1D0) + 
     &          kp4*(-(kp4*p1p3*(p1p3 + p2p3)*2D0) + 
     &             p3p4*(ame**2*(p3p4 - ampi**2*1D0) + 
     &                p1p2*(p3p4 - ampi**2*1D0) + 
     &                p1p3*(p1p3 + kp1*2D0) + 
     &                p2p3*(kp1 - p1p4*1D0 + p1p3*2D0))) + 
     &          kp3*(ame**2*ampi**4 + ampi**4*p1p2 + 
     &             ampi**2*kp1*p2p3 + ampi**2*p1p4*p2p3 - 
     &             ampi**2*p1p3**2*1D0 - 
     &             (ame**2 + p1p2)*p3p4*1D0*
     &              (ampi**2 - kp4*1D0) + 
     &             ampi**2*kp1*p1p3*2D0 + 
     &             ampi**2*p1p3*p1p4*2D0 - 
     &             ampi**2*p1p3*p2p3*2D0 - 
     &             kp4*p1p4*p2p3*2D0 + 
     &             ampi**2*p1p3*p2p4*2D0 - 
     &             kp4*p1p3*p1p4*4D0 - kp4*p1p3*p2p4*4D0)) + 
     &       kp3*(ampi**2*p2p4*
     &           ((kp1 + p1p4)*p2p3 + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             p1p2*(ampi**2 - p3p4*1D0)) + 
     &          kp4**2*p2p3*(p1p2 - ame**2*1D0)*2D0 - 
     &          p1p3**2*p2p4*1D0*(ampi**2 - kp4*2D0) - 
     &          kp4*p2p4*1D0*
     &           (ampi**2*(ame**2 + p1p2)*2D0 + 
     &             p2p3*2D0*(kp1 + p1p4*2D0) - 
     &             (ame**2 + p1p2)*p3p4*3D0) - 
     &          p1p3*1D0*
     &           (ampi**2*
     &              ((kp1 + p1p4)*p2p3 - p2p4**2*1D0 + 
     &                ame**2*(ampi**2 - p3p4*1D0) + 
     &                p1p2*(ampi**2 - p3p4*1D0) - 
     &                p2p4*(kp1 + p1p4 - p2p3*1D0)*2D0) + 
     &             kp4*(-(ampi**2*(ame**2 + p1p2)*2D0) - 
     &                p1p4*p2p3*2D0 + 
     &                (ame**2 + p1p2)*p3p4*3D0 + 
     &                p2p4**2*4D0 + 
     &                p2p4*(kp1 - p2p3*1D0 + p1p4*2D0)*4D0) + 
     &             kp4**2*(-(x12*1D0) + p1p2*4D0 + ame**2*8D0)
     &             ))))/(kp3*kp4**2*(ame**2 + p1p2))
