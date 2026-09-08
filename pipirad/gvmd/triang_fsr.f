        elmat2 = (ame**2*B0ax12m1m2*ep3*2D0*
     &     (ampi**2*kp3**2*p1p3 + ampi**2*kp4**2*p1p3 + 
     &       ampi**2*kp3**2*p2p4 + ampi**2*kp4**2*p2p4 - 
     &       ampi**2*kp3**2*p1p4*1D0 - 
     &       ampi**2*kp4**2*p1p4*1D0 - 
     &       ampi**2*kp3**2*p2p3*1D0 - 
     &       ampi**2*kp4**2*p2p3*1D0 + 
     &       ampi**2*kp1*(kp4**2 - kp3**2*1D0) + 
     &       ampi**2*kp2*(kp3**2 - kp4**2*1D0) - 
     &       kp3*kp4**2*p1p3*2D0 + kp3**2*kp4*p1p4*2D0 + 
     &       kp3*kp4**2*p2p3*2D0 - kp3**2*kp4*p2p4*2D0 + 
     &       kp3*kp4*p3p4*(p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*
     &        2D0))/(kp3**2*kp4**2*(ame**4 - p1p2**2*1D0)) + 
     &  (ame**2*C0aeex12m1em2*ep3*
     &     (ampi**2*kp3**2*p1p4 + ampi**2*kp4**2*p1p4 + 
     &       ampi**2*kp3**2*p2p3 + ampi**2*kp4**2*p2p3 - 
     &       ampi**2*kp3**2*p1p3*1D0 - 
     &       ampi**2*kp4**2*p1p3*1D0 - 
     &       ampi**2*kp3**2*p2p4*1D0 - 
     &       ampi**2*kp4**2*p2p4*1D0 + 
     &       ampi**2*kp2*(kp4**2 - kp3**2*1D0) + 
     &       ampi**2*kp1*(kp3**2 - kp4**2*1D0) + 
     &       kp3*kp4**2*p1p3*2D0 - kp3**2*kp4*p1p4*2D0 - 
     &       kp3*kp4**2*p2p3*2D0 + kp3**2*kp4*p2p4*2D0 + 
     &       kp3*kp4*p3p4*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &        2D0)*(m12 + m22 - p1p2*4D0))/
     &   (kp3**2*kp4**2*(ame**4 - p1p2**2*1D0)) + 
     &  (B0ax15em1*ep2*5.D-1*
     &     (kp1**2*2D0*(kp4*
     &           (p2p3*(kp2 - p1p2*1D0) + 
     &             ame**2*(ampi**2 - p3p4*1D0) + 
     &             ame**2*p1p3*2D0) + 
     &          kp3*(ame**2*ampi**2 + p2p4*(kp2 - p1p2*1D0) - 
     &             ame**2*1D0*(p3p4 - p1p4*2D0 + kp4*4D0))) + 
     &       ame**2*(ame**2 + p1p2 - kp2*1D0)*
     &        (-(kp3**2*(p1p4 + p2p4)*1D0) - 
     &          kp4*1D0*(kp4*(p1p3 + p2p3) + 
     &             p1p3*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &              2D0) + 
     &          kp2*(kp4*(p3p4 - ampi**2*1D0 + p1p3*2D0) + 
     &             kp3*(p3p4 - ampi**2*1D0 + p1p4*2D0 + 
     &                kp4*4D0)) + 
     &          kp3*(p1p4*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &              2D0 + 
     &             kp4*(p1p3 + p1p4 + p2p3 + p2p4 + 
     &                ame**2*4D0 - p1p2*4D0))) - 
     &       kp1*1D0*(-(ame**2*kp3**2*(p1p4 + p2p4)*2D0) + 
     &          kp3*(ame**4*ampi**2 - 
     &             ame**2*ampi**2*p1p2*1D0 - 
     &             ampi**2*p1p2**2*2D0 + ame**4*p1p4*2D0 - 
     &             kp2**2*p1p4*2D0 + ame**2*p1p2*p1p4*2D0 - 
     &             p1p2*p1p4*p2p3*2D0 - p1p2*p1p3*p2p4*2D0 - 
     &             (ame**2 + p1p2)*p3p4*1D0*
     &              (ame**2 - p1p2*2D0) + 
     &             kp4*2D0*
     &              (p1p2**2*2D0 + 
     &                ame**2*
     &                 (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                   ame**2*2D0)) + 
     &             kp2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &                p1p3*p2p4*2D0 + 
     &                p1p2*2D0*
     &                 (ampi**2 + p1p4 - p3p4*1D0 - kp4*2D0)
     &                 + p1p4*2D0*(ame**2 + p2p3 - p2p4*2D0))
     &              + ame**2*p1p3*p1p4*4D0 - 
     &             ame**2*p1p4**2*4D0 - 
     &             ame**2*p1p4*p2p3*4D0 + 
     &             ame**2*p1p4*p2p4*4D0 + p1p2*p1p4*p2p4*4D0)
     &           + kp4*(-(kp2**2*p1p3*2D0) + 
     &             p1p2**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &             kp2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &                p1p4*p2p3*2D0 + 
     &                p1p2*(ampi**2 + p1p3 - p3p4*1D0)*2D0 + 
     &                p1p3*2D0*(ame**2 + p2p4 - p2p3*2D0)) + 
     &             p1p2*(ame**2*(p3p4 - ampi**2*1D0) - 
     &                p1p4*p2p3*2D0 + 
     &                p1p3*2D0*(ame**2 - p2p4*1D0 + p2p3*2D0))
     &               + ame**2*
     &              (ame**2*(ampi**2 - p3p4*1D0) - 
     &                kp4*p2p3*2D0 + 
     &                p1p3*2D0*
     &                 (ame**2 - kp4*1D0 + p1p4*2D0 + 
     &                   p2p3*2D0 - p2p4*2D0) - p1p3**2*4D0)))
     &       ))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)*
     &     (-(ame**4*1D0) + (kp2 - p1p2*1D0)**2 + 
     &       ame**2*kp1*2D0)) + 
     &  (B0ax25em2*ep2*5.D-1*
     &     (-(kp2*1D0*(-(ame**2*kp3**2*(p1p4 + p2p4)*2D0) + 
     &            kp3*(ame**4*ampi**2 - ame**4*p3p4*1D0 + 
     &               p1p2**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &               ame**2*p2p4*2D0*
     &                (ame**2 - p1p3*2D0 + p1p4*2D0 + 
     &                  p2p3*2D0 - p2p4*2D0) + 
     &               kp4*2D0*
     &                (p1p2**2*2D0 + 
     &                  ame**2*
     &                   (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     ame**2*2D0)) + 
     &               p1p2*
     &                (ame**2*(p3p4 - ampi**2*1D0) + 
     &                  p2p4*(ame**2 - p1p3*1D0)*2D0 - 
     &                  p1p4*2D0*(p2p3 - p2p4*2D0))) + 
     &            kp4*(p1p2**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &               p1p2*
     &                (ame**2*(p3p4 - ampi**2*1D0) - 
     &                  p1p3*p2p4*2D0 + 
     &                  p2p3*2D0*
     &                   (ame**2 - p1p4*1D0 + p1p3*2D0)) + 
     &               ame**2*
     &                (ame**2*ampi**2 - ame**2*p3p4*1D0 - 
     &                  kp4*(p1p3 + p2p3)*2D0 + 
     &                  p2p3*2D0*
     &                   (ame**2 + p1p3*2D0 - p1p4*2D0 - 
     &                     p2p3*2D0 + p2p4*2D0))))) + 
     &       kp1**2*(kp4*
     &           (ame**2*(ampi**2 - p3p4*1D0) + 
     &             p2p3*(kp2 - ame**2*1D0)*2D0) + 
     &          kp3*(ame**2*(ampi**2 - p3p4*1D0) + 
     &             p2p4*(kp2 - ame**2*1D0)*2D0 - 
     &             ame**2*kp4*4D0)) + 
     &       kp2**2*2D0*(kp4*
     &           (-(p1p2*p1p3*1D0) + 
     &             ame**2*(ampi**2 - p3p4*1D0 + p2p3*2D0)) + 
     &          kp3*(-(p1p2*p1p4*1D0) + 
     &             ame**2*(ampi**2 - p3p4*1D0 + p2p4*2D0) - 
     &             ame**2*kp4*4D0)) + 
     &       ame**2*(ame**2 + p1p2)*
     &        (-(kp3**2*(p1p4 + p2p4)*1D0) - 
     &          kp4*1D0*(kp4*(p1p3 + p2p3) + 
     &             p2p3*(p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*
     &              2D0) + 
     &          kp3*(p2p4*(p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*
     &              2D0 + 
     &             kp4*(p1p3 + p1p4 + p2p3 + p2p4 + 
     &                ame**2*4D0 - p1p2*4D0))) + 
     &       kp1*(ame**2*kp3**2*(p1p4 + p2p4) + 
     &          kp4*(kp2**2*p1p3*2D0 + 
     &             kp2*(ame**2*(ampi**2 - p3p4*1D0) - 
     &                p1p3*p2p4*2D0 - 
     &                p1p2*(ampi**2 + p2p3 - p3p4*1D0)*2D0 - 
     &                p2p3*2D0*(ame**2 + p1p4 - p1p3*2D0)) + 
     &             ame**2*
     &              (kp4*(p1p3 + p2p3) + ame**2*p3p4 - 
     &                ame**2*ampi**2*1D0 + 
     &                p2p3*
     &                 (ame**2 + p1p4 + p2p3 - p1p3*1D0 - 
     &                   p2p4*1D0)*2D0 + 
     &                p1p2*(p3p4 - ampi**2*1D0 + p2p3*2D0)))
     &           - kp3*1D0*
     &           (-(kp2**2*p1p4*2D0) + 
     &             kp2*(ame**2*(p3p4 - ampi**2*1D0) + 
     &                (ame**2 + p1p3)*p2p4*2D0 + 
     &                p1p2*2D0*
     &                 (ampi**2 + p2p4 - p3p4*1D0 - kp4*2D0)
     &                 + p1p4*2D0*(p2p3 - p2p4*2D0)) + 
     &             ame**2*
     &              (ame**2*ampi**2 + 
     &                kp4*(p1p3 + p1p4 + p2p3 + p2p4) - 
     &                ame**2*p3p4*1D0 - 
     &                p2p4*
     &                 (ame**2 + p1p3 + p2p4 - p1p4*1D0 - 
     &                   p2p3*1D0)*2D0 + 
     &                p1p2*
     &                 (ampi**2 - p3p4*1D0 - p2p4*2D0 - 
     &                   kp4*8D0))))))/
     &   (kp2*kp3*kp4*(ame**2 + p1p2)*
     &     (-(ame**4*1D0) + (kp1 - p1p2*1D0)**2 + 
     &       ame**2*kp2*2D0)) + 
     &  (ame**2*B0aeem1*5.D-1*
     &     (ampi**2*ep3*kp1**3*kp2*(kp3**2 - kp4**2*1D0)*
     &        2D0 + ame**2*ampi**2*ep3*kp2**3*
     &        (kp4**2 - kp3**2*1D0)*4D0 + 
     &       kp1**2*(ampi**2*ep3*kp2**2*(kp4**2 - kp3**2*1D0)*
     &           2D0 + ep3*kp2*2D0*
     &           (kp3*kp4*
     &              (kp4*(p1p3 - p2p3*1D0) + 
     &                p3p4*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)
     &                )*2D0 + 
     &             ampi**2*kp4**2*
     &              (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0 + 
     &                p1p2*2D0) + 
     &             kp3**2*
     &              (ampi**2*p1p4 + ampi**2*p2p3 - 
     &                ampi**2*p1p3*1D0 - ampi**2*p1p2*2D0 - 
     &                kp4*p1p4*2D0 - 
     &                p2p4*1D0*(ampi**2 - kp4*2D0))) - 
     &          ep2*kp3*kp4*1D0*(ame**2 - p1p2*1D0)*
     &           (kp4*(p3p4 - ampi**2*1D0 + p2p3*2D0) + 
     &             kp3*(p3p4 - ampi**2*1D0 + p2p4*2D0 + 
     &                kp4*4D0))) + 
     &       ep2*kp3*kp4*(ame**4 - p1p2**2*1D0)*
     &        (-(kp3**2*(p1p4 + p2p4)*1D0) - 
     &          kp4*1D0*(kp4*(p1p3 + p2p3) + 
     &             p2p3*(p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*
     &              2D0) + 
     &          kp3*(p2p4*(p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*
     &              2D0 + 
     &             kp4*(p1p3 + p1p4 + p2p3 + p2p4 + 
     &                ame**2*4D0 - p1p2*4D0))) + 
     &       kp2**2*2D0*(ampi**2*ep3*kp4**2*
     &           (p1p2**2 - 
     &             ame**2*1D0*
     &              (ame**2 + p1p3*2D0 - p1p4*2D0 - 
     &                p2p3*2D0 + p2p4*2D0)) + 
     &          kp3*kp4*(ame**2*ep3*p3p4*
     &              (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*4D0 + 
     &             kp4*(-(ame**2*ep3*p2p3*4D0) + 
     &                p1p3*
     &                 (ep2*p1p2 - ame**2*1D0*(ep2 - ep3*4D0))
     &                )) + 
     &          kp3**2*(ep2*kp4*p1p2*p1p4 - 
     &             ampi**2*ep3*p1p2**2*1D0 + 
     &             ame**2*
     &              (-(ampi**2*ep3*p1p3*2D0) - 
     &                ep3*p2p4*2D0*(ampi**2 - kp4*2D0) + 
     &                ampi**2*ep3*(ame**2 + p2p3*2D0) + 
     &                p1p4*
     &                 (ampi**2*ep3*2D0 - 
     &                   kp4*1D0*(ep2 + ep3*4D0))))) + 
     &       kp2*(ame**2 - p1p2*1D0)*
     &        (ampi**2*ep3*kp4**2*(ame**2 + p1p2)*
     &           (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*2D0 + 
     &          kp3*kp4*(-(ep3*(ame**2 + p1p2)*p3p4*
     &                (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*4D0)
     &               + kp4*
     &              (ame**2*ampi**2*ep2 - 
     &                ame**2*ep2*p3p4*1D0 + 
     &                p2p3*2D0*
     &                 (ep2*p1p4 - ame**2*1D0*(ep2 - ep3*2D0))
     &                  + p1p2*
     &                 (ampi**2*ep2 - ep2*p3p4*1D0 - 
     &                   p2p3*2D0*(ep2 - ep3*2D0) - 
     &                   ep3*p1p3*4D0) + 
     &                p1p3*
     &                 (ep2*p2p4*2D0 - ame**2*ep3*4D0 - 
     &                   ep2*p2p3*4D0))) + 
     &          kp3**2*(ampi**2*ep3*(ame**2 + p1p2)*
     &              (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*2D0 + 
     &             kp4*(ame**2*ampi**2*ep2 - 
     &                ame**2*ep2*p3p4*1D0 + 
     &                p1p4*2D0*
     &                 (ep2*p2p3 + ame**2*ep3*2D0 - 
     &                   ep2*p2p4*2D0) - 
     &                p2p4*2D0*
     &                 (-(ep2*p1p3*1D0) + 
     &                   ame**2*(ep2 + ep3*2D0)) + 
     &                p1p2*
     &                 (ampi**2*ep2 - ep2*p3p4*1D0 - 
     &                   p2p4*2D0*(ep2 + ep3*2D0) + 
     &                   ep3*p1p4*4D0)) - 
     &             ame**2*ep2*kp4**2*8D0)) + 
     &       kp1*(ampi**2*ep3*kp2**2*(ame**2 + p1p2)*
     &           (kp3**2 - kp4**2*1D0)*4D0 + 
     &          ep2*kp3*kp4*(ame**2 - p1p2*1D0)*
     &           (kp3**2*(p1p4 + p2p4) + 
     &             kp4*(kp4*(p1p3 + p2p3) + ame**2*p3p4 - 
     &                ame**2*ampi**2*1D0 + 
     &                p2p3*
     &                 (ame**2 + p1p4 + p2p3 - p1p3*1D0 - 
     &                   p2p4*1D0)*2D0 + 
     &                p1p2*(p3p4 - ampi**2*1D0 + p2p3*2D0)) - 
     &             kp3*1D0*
     &              (ame**2*ampi**2 + 
     &                kp4*(p1p3 + p1p4 + p2p3 + p2p4) - 
     &                ame**2*p3p4*1D0 - 
     &                p2p4*
     &                 (ame**2 + p1p3 + p2p4 - p1p4*1D0 - 
     &                   p2p3*1D0)*2D0 + 
     &                p1p2*
     &                 (ampi**2 - p3p4*1D0 - p2p4*2D0 - 
     &                   kp4*8D0))) + 
     &          kp2*(ampi**2*ep3*kp4**2*2D0*
     &              (ame**4 - p1p2**2*1D0 + 
     &                p1p2*
     &                 (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*2D0
     &                ) + 
     &             kp3*kp4*
     &              (ep3*p1p2*p3p4*
     &                 (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*8D0
     &                  + kp4*
     &                 (ame**2*ep2*(ampi**2 - p3p4*1D0) + 
     &                   p1p2*
     &                    (ep2*p3p4 - ampi**2*ep2*1D0 - 
     &                     ep3*p1p3*8D0 + ep3*p2p3*8D0))) + 
     &             kp3**2*
     &              (ampi**2*ep3*2D0*
     &                 (p1p2**2 - ame**4*1D0 + 
     &                   p1p2*
     &                    (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &                    2D0) + 
     &                ep2*kp4**2*(p1p2 - ame**2*1D0)*4D0 + 
     &                kp4*
     &                 (ame**2*ep2*(ampi**2 - p3p4*1D0) + 
     &                   p1p2*
     &                    (ep2*p3p4 - ampi**2*ep2*1D0 + 
     &                     ep3*p1p4*8D0 - ep3*p2p4*8D0)))))))/
     &   (kp2*kp3**2*kp4**2*(ame**4 - p1p2**2*1D0)*
     &     (-(ame**4*1D0) + (kp1 - p1p2*1D0)**2 + 
     &       ame**2*kp2*2D0)) + 
     &  (ame**2*B0aeem2*5.D-1*
     &     (ame**2*ampi**2*ep3*kp1**3*(kp3**2 - kp4**2*1D0)*
     &        4D0 + ep2*kp3*kp4*(ame**2 + p1p2 - kp2*1D0)*
     &        (ame**2 - p1p2*1D0)*
     &        (-(kp3**2*(p1p4 + p2p4)*1D0) - 
     &          kp4*1D0*(kp4*(p1p3 + p2p3) + 
     &             p1p3*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &              2D0) + 
     &          kp2*(kp4*(p3p4 - ampi**2*1D0 + p1p3*2D0) + 
     &             kp3*(p3p4 - ampi**2*1D0 + p1p4*2D0 + 
     &                kp4*4D0)) + 
     &          kp3*(p1p4*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &              2D0 + 
     &             kp4*(p1p3 + p1p4 + p2p3 + p2p4 + 
     &                ame**2*4D0 - p1p2*4D0))) - 
     &       kp1**2*2D0*(ampi**2*ep3*kp4**2*
     &           (kp2**2 + p1p2**2 - ame**4*1D0 - 
     &             kp2*(ame**2 + p1p2)*2D0 + 
     &             ame**2*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &              2D0) + 
     &          kp3**2*(ame**4*ampi**2*ep3 - 
     &             ampi**2*ep3*kp2**2*1D0 + 
     &             ampi**2*ep3*kp2*(ame**2 + p1p2)*2D0 - 
     &             ep3*1D0*
     &              (ampi**2*p1p2**2 + 
     &                ame**2*2D0*
     &                 (ampi**2*p2p3 - ampi**2*p1p3*1D0 + 
     &                   p1p4*(ampi**2 - kp4*2D0))) + 
     &             p2p4*(ame**2*ampi**2*ep3*2D0 + 
     &                kp4*
     &                 (-(ep2*p1p2*1D0) + 
     &                   ame**2*(ep2 - ep3*4D0)))) + 
     &          kp3*kp4*(ame**2*ep3*p3p4*
     &              (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*4D0 + 
     &             kp4*(-(ame**2*ep3*p1p3*4D0) + 
     &                p2p3*
     &                 (-(ep2*p1p2*1D0) + 
     &                   ame**2*(ep2 + ep3*4D0))))) + 
     &       kp1*(ampi**2*ep3*kp2**3*(kp4**2 - kp3**2*1D0)*
     &           2D0 + ep3*kp2**2*2D0*
     &           (kp3*kp4*
     &              (kp4*(p1p3 - p2p3*1D0) + 
     &                p3p4*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)
     &                )*2D0 + 
     &             ampi**2*kp4**2*
     &              (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0 - 
     &                p1p2*2D0) + 
     &             kp3**2*
     &              (ampi**2*p1p4 + ampi**2*p2p3 - 
     &                ampi**2*p1p3*1D0 + ampi**2*p1p2*2D0 - 
     &                kp4*p1p4*2D0 - 
     &                p2p4*1D0*(ampi**2 - kp4*2D0))) + 
     &          (ame**2 - p1p2*1D0)*
     &           (ampi**2*ep3*kp4**2*(ame**2 + p1p2)*
     &              (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*2D0 + 
     &             kp3*kp4*
     &              (-(ep3*(ame**2 + p1p2)*p3p4*
     &                   (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &                   4D0) + 
     &                kp4*
     &                 (ame**2*ampi**2*ep2 - 
     &                   ame**2*ep2*p3p4*1D0 + 
     &                   p2p3*2D0*
     &                    (ep2*p1p4 + ame**2*ep3*2D0) - 
     &                   p1p3*2D0*
     &                    (-(ep2*p2p4*1D0) + ep2*p2p3*2D0 + 
     &                     ame**2*(ep2 + ep3*2D0)) + 
     &                   p1p2*
     &                    (ampi**2*ep2 - ep2*p3p4*1D0 - 
     &                     p1p3*2D0*(ep2 + ep3*2D0) + 
     &                     ep3*p2p3*4D0))) + 
     &             kp3**2*
     &              (ampi**2*ep3*(ame**2 + p1p2)*
     &                 (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*2D0
     &                  + kp4*
     &                 (ame**2*ampi**2*ep2 - 
     &                   ame**2*ep2*p3p4*1D0 + 
     &                   p2p4*2D0*
     &                    (ep2*p1p3 - ame**2*ep3*2D0) + 
     &                   p1p4*
     &                    (-(ame**2*2D0*(ep2 - ep3*2D0)) + 
     &                     ep2*2D0*(p2p3 - p2p4*2D0)) + 
     &                   p1p2*
     &                    (ampi**2*ep2 - ep2*p3p4*1D0 - 
     &                     p1p4*2D0*(ep2 - ep3*2D0) - 
     &                     ep3*p2p4*4D0)) - 
     &                ame**2*ep2*kp4**2*8D0)) + 
     &          kp2*(ampi**2*ep3*kp4**2*2D0*
     &              (p1p2**2 - ame**4*1D0 + 
     &                p1p2*
     &                 (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*2D0
     &                ) + 
     &             kp3*kp4*
     &              (ep3*p1p2*p3p4*
     &                 (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*8D0
     &                  + kp4*
     &                 (ame**2*ep2*(ampi**2 - p3p4*1D0) + 
     &                   p1p2*
     &                    (ep2*p3p4 - ampi**2*ep2*1D0 - 
     &                     ep3*p1p3*8D0 + ep3*p2p3*8D0))) + 
     &             kp3**2*
     &              (ampi**2*ep3*2D0*
     &                 (ame**4 - p1p2**2*1D0 + 
     &                   p1p2*
     &                    (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &                    2D0) + 
     &                ep2*kp4**2*(p1p2 - ame**2*1D0)*4D0 + 
     &                kp4*
     &                 (ame**2*ep2*(ampi**2 - p3p4*1D0) + 
     &                   p1p2*
     &                    (ep2*p3p4 - ampi**2*ep2*1D0 + 
     &                     ep3*p1p4*8D0 - ep3*p2p4*8D0)))))))/
     &   (kp1*kp3**2*kp4**2*(ame**4 - p1p2**2*1D0)*
     &     (-(ame**4*1D0) + (kp2 - p1p2*1D0)**2 + 
     &       ame**2*kp1*2D0)) - 
     &  (C0aex15x34m2em1*ep2*5.D-1*
     &     (ame**2*kp1**3*(kp4*p2p3 + kp3*p2p4)*4D0 - 
     &       kp1**2*2D0*(kp4*
     &           (kp2*(ame**2*(ampi**2 - p3p4*1D0) + 
     &                p2p3*(m22 - ame**2*2D0)) + 
     &             p1p2*(ame**2*(ampi**2 - p3p4*1D0) - 
     &                ame**2*p1p3*2D0 + 
     &                p2p3*(-(m22*1D0) + ame**2*2D0)) + 
     &             ame**2*
     &              (p2p3*
     &                 (-(m12*1D0) + ame**2*2D0 + p1p4*2D0) + 
     &                p1p3*2D0*
     &                 (ame**2 + m22 + p2p4 - p2p3*2D0) + 
     &                (ampi**2 - p3p4*1D0)*(m22 + ame**2*3D0))
     &             ) + kp3*
     &           (ame**2*ampi**2*p1p2 - ame**2*m12*p2p4*1D0 - 
     &             m22*p1p2*p2p4*1D0 + ame**4*p1p4*2D0 + 
     &             ame**2*m22*p1p4*2D0 - 
     &             ame**2*p1p2*p1p4*2D0 + 
     &             ame**2*p1p4*p2p3*2D0 + ame**4*p2p4*2D0 + 
     &             ame**2*p1p2*p2p4*2D0 + 
     &             ame**2*p1p3*p2p4*2D0 + 
     &             ame**2*ampi**2*(m22 + ame**2*3D0) - 
     &             ame**2*p3p4*1D0*
     &              (m22 + p1p2 + ame**2*3D0) - 
     &             ame**2*p1p4*p2p4*4D0 - 
     &             ame**2*kp4*4D0*(m22 + ame**2*4D0) + 
     &             kp2*(ame**2*(ampi**2 - p3p4*1D0) + 
     &                p2p4*(m22 - ame**2*2D0) - ame**2*kp4*4D0
     &                ))) - 
     &       ame**2*1D0*(ame**2 + p1p2 - kp2*1D0)*
     &        (m12 + m22 + kp2*4D0 - p1p2*4D0)*
     &        (-(kp3**2*(p1p4 + p2p4)*1D0) - 
     &          kp4*1D0*(kp4*(p1p3 + p2p3) + 
     &             p1p3*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &              2D0) + 
     &          kp2*(kp4*(p3p4 - ampi**2*1D0 + p1p3*2D0) + 
     &             kp3*(p3p4 - ampi**2*1D0 + p1p4*2D0 + 
     &                kp4*4D0)) + 
     &          kp3*(p1p4*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &              2D0 + 
     &             kp4*(p1p3 + p1p4 + p2p3 + p2p4 + 
     &                ame**2*4D0 - p1p2*4D0))) + 
     &       kp1*(-(ame**2*kp3**2*(p1p4 + p2p4)*
     &             (ame**2 + kp2 + m22 - p1p2*1D0)*2D0) + 
     &          kp4*(-(kp2**2*2D0*
     &                (m22*p1p3 + 
     &                  ame**2*(ampi**2 - p3p4*1D0)*3D0)) - 
     &             p1p2**2*2D0*
     &              (m22*(ampi**2 - p3p4*1D0) + 
     &                ame**2*p1p3*4D0) + 
     &             ame**2*
     &              (-(p2p3*2D0*
     &                   (kp4*(ame**2 + m22) + 
     &                     p1p4*(m12 - ame**2*2D0))) + 
     &                p1p3*2D0*
     &                 (ame**2*(m12 + m22) - 
     &                   kp4*(ame**2 + m22)*1D0 + 
     &                   (ame**2 + m22)*p1p4*2D0 + 
     &                   p2p3*(m12 + m22 - ame**2*1D0)*2D0 - 
     &                   p2p4*1D0*(m12 + m22*2D0)) - 
     &                (ame**2 + m22)*p1p3**2*4D0 + 
     &                ame**2*p3p4*
     &                 (m12 - m22*1D0 - ame**2*4D0) + 
     &                ame**2*ampi**2*
     &                 (m22 - m12*1D0 + ame**2*4D0)) + 
     &             kp2*(ame**2*p3p4*
     &                 (m12 + m22 + ame**2*2D0) - 
     &                ame**2*ampi**2*1D0*
     &                 (m12 + m22 + ame**2*2D0) + 
     &                p2p3*2D0*
     &                 (-(ame**2*kp4*1D0) + 
     &                   p1p4*(m22 - ame**2*2D0)) - 
     &                ame**2*p1p3**2*4D0 + 
     &                p1p3*2D0*
     &                 (-(ame**2*kp4*1D0) + ame**2*p1p4*2D0 + 
     &                   p2p3*2D0*(-(m22*1D0) + ame**2*3D0) + 
     &                   p2p4*(m22 - ame**2*4D0) + 
     &                   ame**2*(m22 + ame**2*4D0)) + 
     &                p1p2*2D0*
     &                 ((ampi**2 - p3p4*1D0)*
     &                    (m22 + ame**2*3D0) + 
     &                   p1p3*(m22 + ame**2*4D0))) + 
     &             p1p2*(p2p3*2D0*
     &                 (ame**2*kp4 + 
     &                   p1p4*(-(m22*1D0) + ame**2*2D0)) + 
     &                ame**2*p1p3**2*4D0 + 
     &                ame**2*p3p4*(m12 + m22 - ame**2*4D0) + 
     &                ame**2*ampi**2*
     &                 (-(m12*1D0) - m22*1D0 + ame**2*4D0) + 
     &                p1p3*2D0*
     &                 (ame**2*kp4 - ame**2*p1p4*2D0 + 
     &                   p2p3*2D0*(m22 - ame**2*3D0) + 
     &                   ame**2*(m12 + m22 - ame**2*4D0) + 
     &                   p2p4*(-(m22*1D0) + ame**2*4D0)))) + 
     &          kp3*(-(ame**2*ampi**2*m12*p1p2*1D0) - 
     &             ame**2*ampi**2*m22*p1p2*1D0 - 
     &             ampi**2*m22*p1p2**2*2D0 + 
     &             ame**4*m12*p1p4*2D0 + 
     &             ame**4*m22*p1p4*2D0 + 
     &             ame**2*m12*p1p2*p1p4*2D0 + 
     &             ame**2*m22*p1p2*p1p4*2D0 - 
     &             ame**2*m12*p1p4*p2p3*2D0 - 
     &             m22*p1p2*p1p4*p2p3*2D0 - 
     &             ame**2*m12*p1p3*p2p4*2D0 - 
     &             m22*p1p2*p1p3*p2p4*2D0 + 
     &             ame**4*ampi**2*p1p2*4D0 + 
     &             ame**4*p1p3*p1p4*4D0 + 
     &             ame**2*m22*p1p3*p1p4*4D0 - 
     &             ame**2*p1p2*p1p3*p1p4*4D0 - 
     &             ame**4*p1p4**2*4D0 - 
     &             ame**2*m22*p1p4**2*4D0 + 
     &             ame**2*p1p2*p1p4**2*4D0 - 
     &             ame**2*m22*p1p4*p2p3*4D0 + 
     &             ame**4*p1p3*p2p4*4D0 + 
     &             ame**2*p1p2*p1p3*p2p4*4D0 - 
     &             ame**4*p1p4*p2p4*4D0 + 
     &             ame**2*m12*p1p4*p2p4*4D0 + 
     &             ame**2*m22*p1p4*p2p4*4D0 + 
     &             m22*p1p2*p1p4*p2p4*4D0 + 
     &             ame**4*ampi**2*
     &              (m22 - m12*1D0 + ame**2*4D0) - 
     &             (ame**2 + p1p2)*p3p4*1D0*
     &              (-(m22*p1p2*2D0) + 
     &                ame**2*(m22 - m12*1D0 + ame**2*4D0)) - 
     &             ame**4*p1p2*p1p4*8D0 - 
     &             ame**2*p1p2**2*p1p4*8D0 + 
     &             ame**2*p1p2*p1p4*p2p3*8D0 - 
     &             ame**2*p1p2*p1p4*p2p4*1.2D1 + 
     &             kp2*(ame**2*p3p4*
     &                 (m12 + m22 + ame**2*2D0) - 
     &                ame**2*ampi**2*1D0*
     &                 (m12 + m22 + ame**2*2D0) + 
     &                p2p4*2D0*
     &                 (p1p3*(m22 - ame**2*2D0) + 
     &                   p1p4*2D0*(-(m22*1D0) + ame**2*3D0))
     &                 + p1p2*2D0*
     &                 ((ampi**2 - p3p4*1D0)*
     &                    (m22 + ame**2*3D0) + 
     &                   p1p4*(m22 + ame**2*4D0)) + 
     &                p1p4*2D0*
     &                 (p2p3*(m22 - ame**2*4D0) + 
     &                   ame**2*
     &                    (m22 + p1p3*2D0 - p1p4*2D0 + 
     &                     ame**2*4D0)) + 
     &                kp4*2D0*
     &                 (-(p1p2*2D0*(m22 + ame**2*1.D1)) + 
     &                   ame**2*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     m12*2D0 + ame**2*1.2D1))) + 
     &             kp4*2D0*
     &              (ame**2*(ame**2 + m22)*
     &                 (p1p3 + p1p4 + p2p3 + p2p4) + 
     &                ame**4*2D0*
     &                 (m22 + m12*2D0 - ame**2*4D0) + 
     &                p1p2**2*2D0*(m22 + ame**2*4D0) - 
     &                ame**2*p1p2*1D0*
     &                 (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                   ame**2*1.6D1)) + 
     &             kp2**2*
     &              (-(m22*p1p4*2D0) + 
     &                ame**2*(p3p4 - ampi**2*1D0)*6D0 + 
     &                ame**2*kp4*2.4D1)))))/
     &   (kp1*kp3*kp4*(ame**2 + p1p2)*
     &     (-(ame**4*1D0) + (kp2 - p1p2*1D0)**2 + 
     &       ame**2*kp1*2D0)) + 
     &  (B0ax34m1m2*ep2*(-(kp1**4*1D0*
     &          (kp4*(kp2**2*p2p3 + 
     &               kp2*
     &                (p2p3*(ame**2 - p1p2*1D0) + 
     &                  ame**2*(ampi**2 - p3p4*1D0) + 
     &                  ame**2*p1p3*2D0) + 
     &               ame**4*2D0*
     &                (ampi**2 - p3p4*1D0 - p2p3*2D0)) + 
     &            kp3*(kp2**2*p2p4 + 
     &               ame**4*2D0*
     &                (ampi**2 - p3p4*1D0 - p2p4*2D0 - 
     &                  kp4*4D0) + 
     &               kp2*
     &                (ame**2*ampi**2 - ame**2*p3p4*1D0 + 
     &                  p2p4*(ame**2 - p1p2*1D0) + 
     &                  ame**2*p1p4*2D0 - ame**2*kp4*4D0))))
     &        + ame**2*kp2*(ame**2 + p1p2 - kp2*1D0)*
     &        (ame**4 - p1p2**2*1D0 - ame**2*kp2*2D0)*
     &        (-(kp3**2*(p1p4 + p2p4)*1D0) - 
     &          kp4*1D0*(kp4*(p1p3 + p2p3) + 
     &             p1p3*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &              2D0) + 
     &          kp2*(kp4*(p3p4 - ampi**2*1D0 + p1p3*2D0) + 
     &             kp3*(p3p4 - ampi**2*1D0 + p1p4*2D0 + 
     &                kp4*4D0)) + 
     &          kp3*(p1p4*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &              2D0 + 
     &             kp4*(p1p3 + p1p4 + p2p3 + p2p4 + 
     &                ame**2*4D0 - p1p2*4D0))) + 
     &       kp1**3*(-(ame**2*kp3**2*(p1p4 + p2p4)*1D0*
     &             (kp2 + ame**2*2D0)) - 
     &          kp4*1D0*(kp2**3*(p1p3 + p2p3) + 
     &             ame**2*p1p2**2*
     &              (ampi**2 - p3p4*1D0 - p2p3*2D0) + 
     &             ame**4*p1p2*2D0*
     &              (p3p4 - ampi**2*1D0 + p2p3*2D0) + 
     &             ame**4*
     &              (kp4*(p1p3 + p2p3)*2D0 + 
     &                ame**2*(p3p4 - ampi**2*1D0 + p2p3*2D0)*
     &                 3D0 + 
     &                p2p3*
     &                 (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*4D0
     &                ) + 
     &             kp2**2*
     &              (ame**2*ampi**2*2D0 - ame**2*p3p4*2D0 - 
     &                p2p3*1D0*(p1p4 + ame**2*2D0) + 
     &                p1p3*(ame**2 - p2p4*1D0 + p2p3*2D0) - 
     &                p1p2*1D0*
     &                 (ampi**2 + p1p3 - p3p4*1D0 + p2p3*4D0))
     &               + kp2*
     &              (p1p2**2*
     &                 (ampi**2 - p3p4*1D0 + p2p3*3D0) + 
     &                ame**2*
     &                 (p1p3**2*2D0 + 
     &                   ame**2*(ampi**2 - p3p4*1D0)*2D0 + 
     &                   p1p3*
     &                    (kp4 + p2p4 - (ame**2 + p1p4)*2D0)
     &                    + p2p3*(kp4 - p1p4*1D0 - ame**2*3D0)
     &                   ) + 
     &                p1p2*
     &                 (p2p3*(p1p4 + ame**2*4D0) + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*5D0 + 
     &                   p1p3*(p2p4 - p2p3*2D0 - ame**2*6D0)))
     &             ) + kp3*
     &           (-(kp2**3*(p1p4 + p2p4)*1D0) + 
     &             kp2**2*
     &              (ampi**2*(p1p2 - ame**2*2D0) + 
     &                p3p4*(-(p1p2*1D0) + ame**2*2D0) + 
     &                p1p4*
     &                 (p1p2 + p2p3 - ame**2*1D0 - p2p4*2D0)
     &                 + p2p4*
     &                 (p1p3 + ame**2*2D0 + p1p2*4D0) + 
     &                kp4*(-(p1p2*2D0) + ame**2*6D0)) + 
     &             ame**2*
     &              (ame**2*p1p2*2D0*
     &                 (ampi**2 - p3p4*1D0 - p2p4*2D0) + 
     &                p1p2**2*
     &                 (p3p4 - ampi**2*1D0 + p2p4*2D0) + 
     &                ame**4*ampi**2*3D0 - ame**4*p3p4*3D0 - 
     &                ame**2*p2p4*2D0*
     &                 (p1p3*2D0 - p1p4*2D0 - p2p3*2D0 + 
     &                   p2p4*2D0 + ame**2*3D0) + 
     &                kp4*2D0*
     &                 (p1p2**2*2D0 + 
     &                   ame**2*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 - 
     &                     ame**2*2D0) - ame**2*p1p2*8D0)) + 
     &             kp2*(p1p2**2*
     &                 (p3p4 - ampi**2*1D0 - p2p4*3D0) + 
     &                ame**2*
     &                 (-(p1p4**2*2D0) + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &                   p1p4*
     &                    (-(p2p3*1D0) + (ame**2 + p1p3)*2D0)
     &                    + p2p4*(p1p3 + ame**2*3D0)) + 
     &                p1p2*
     &                 (-(p2p4*1D0*(p1p3 + ame**2*4D0)) + 
     &                   ame**2*(ampi**2 - p3p4*1D0)*5D0 + 
     &                   p1p4*
     &                    (-(p2p3*1D0) + p2p4*2D0 + 
     &                     ame**2*6D0)) + 
     &                kp4*
     &                 (p1p2**2*2D0 + 
     &                   ame**2*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     ame**2*1.D1) - ame**2*p1p2*2.D1))))
     &         - kp1**2*1D0*
     &        (-(ame**2*kp3**2*(p1p4 + p2p4)*1D0*
     &             (-(p1p2**2*1D0) - 
     &               kp2*1D0*(ame**2 + kp2*2D0) + 
     &               ame**4*3D0 + p1p2*(ame**2*2D0 + kp2*5D0))
     &             ) + kp4*
     &           (kp2**4*p1p3 + 
     &             ame**2*(ame**2 + p1p2)*
     &              (p1p2**2*
     &                 (p3p4 - ampi**2*1D0 + p2p3*2D0) + 
     &                p1p2*
     &                 (kp4*(p1p3 + p2p3) + 
     &                   p2p3*
     &                    (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*
     &                    2D0) + 
     &                ame**2*
     &                 (ame**2*ampi**2 - ame**2*p3p4*1D0 - 
     &                   kp4*(p1p3 + p2p3)*3D0 - 
     &                   p2p3*2D0*
     &                    (ame**2 - p1p3*3D0 + p1p4*3D0 + 
     &                     p2p3*3D0 - p2p4*3D0))) + 
     &             kp2**3*
     &              (p2p3*(ame**2 - p1p4*1D0) + 
     &                ame**2*ampi**2*2D0 - ame**2*p3p4*2D0 - 
     &                p1p3*1D0*
     &                 (p2p4 + ame**2*2D0 - p2p3*2D0) - 
     &                p1p2*1D0*
     &                 (ampi**2 + p2p3 - p3p4*1D0 + p1p3*4D0))
     &               + kp2**2*
     &              (ame**4*ampi**2*2D0 + 
     &                p1p2**2*
     &                 (p1p3*3D0 + p2p3*3D0 + ampi**2*4D0 - 
     &                   p3p4*4D0) + 
     &                p1p2*2D0*
     &                 (p2p3*(ame**2 + p1p4*2D0) + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*3D0 + 
     &                   p1p3*(ame**2 + p2p4*2D0 - p2p3*4D0))
     &                 + ame**2*
     &                 (p1p3**2*2D0 - ame**2*p3p4*2D0 + 
     &                   p2p3*
     &                    (kp4*2D0 + p1p4*2D0 + p2p3*2D0 - 
     &                     p2p4*2D0 + ame**2*3D0) + 
     &                   p1p3*
     &                    (kp4*2D0 - p1p4*2D0 + p2p4*2D0 + 
     &                     ame**2*3D0 - p2p3*4D0))) + 
     &             kp2*(p1p2**3*
     &                 (-(p2p3*2D0) - ampi**2*3D0 + p3p4*3D0)
     &                 + p1p2**2*
     &                 (-(p2p3*(p1p4 + ame**2*2D0)*3D0) + 
     &                   ame**2*(ampi**2 - p3p4*1D0)*4D0 + 
     &                   p1p3*
     &                    (-(p2p4*3D0) + (ame**2 + p2p3)*6D0))
     &                  + ame**4*
     &                 (-(p1p3**2*2D0) + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &                   p2p3**2*4D0 + 
     &                   p2p3*
     &                    (kp4 - ame**2*2D0 - p2p4*4D0 + 
     &                     p1p4*7D0) + 
     &                   p1p3*
     &                    (kp4 + p2p4 - ame**2*2D0 + 
     &                     p1p4*2D0 - p2p3*8D0)) + 
     &                ame**2*p1p2*
     &                 (ame**2*ampi**2*5D0 - 
     &                   ame**2*p3p4*5D0 - p1p3**2*6D0 - 
     &                   p2p3*1D0*
     &                    (p1p4*4D0 + p2p3*4D0 - p2p4*4D0 + 
     &                     kp4*5D0 + ame**2*6D0) + 
     &                   p1p3*
     &                    (ame**2*4D0 - kp4*5D0 + p1p4*6D0 - 
     &                     p2p4*6D0 + p2p3*1.D1)))) + 
     &          kp3*(kp2**4*p1p4 + 
     &             kp2**3*
     &              (p2p4*(ame**2 - p1p2*1D0 - p1p3*1D0) + 
     &                p3p4*(p1p2 - ame**2*2D0) + 
     &                ampi**2*(-(p1p2*1D0) + ame**2*2D0) + 
     &                kp4*2D0*(p1p2 - ame**2*3D0) - 
     &                p1p4*1D0*
     &                 (p2p3 + ame**2*2D0 - p2p4*2D0 + 
     &                   p1p2*4D0)) + 
     &             kp2**2*
     &              (ame**4*ampi**2*2D0 + 
     &                p1p2**2*
     &                 (p1p4*3D0 + p2p4*3D0 + ampi**2*4D0 - 
     &                   p3p4*4D0) + 
     &                p1p2*2D0*
     &                 (p2p4*(ame**2 + p1p3*2D0) + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*3D0 + 
     &                   p1p4*(ame**2 + p2p3*2D0 - p2p4*4D0))
     &                 + ame**2*
     &                 (p1p4**2*2D0 - ame**2*p3p4*2D0 + 
     &                   p2p4*
     &                    (p1p3*2D0 - p2p3*2D0 + p2p4*2D0 + 
     &                     ame**2*3D0) + 
     &                   p1p4*
     &                    (-(p1p3*2D0) + p2p3*2D0 + 
     &                     ame**2*3D0 - p2p4*4D0)) - 
     &                kp4*2D0*
     &                 (p1p2**2*4D0 + 
     &                   ame**2*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     ame**2*8D0) - ame**2*p1p2*1.2D1))
     &              + ame**2*(ame**2 + p1p2)*
     &              (ame**4*ampi**2 - ame**4*p3p4*1D0 + 
     &                p1p2*p2p4*
     &                 (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*2D0
     &                  + p1p2**2*
     &                 (p3p4 - ampi**2*1D0 + p2p4*2D0) - 
     &                ame**2*p2p4*2D0*
     &                 (ame**2 + p1p3*3D0 - p1p4*3D0 - 
     &                   p2p3*3D0 + p2p4*3D0) + 
     &                kp4*
     &                 (ame**2*(p1p3 + p1p4 + p2p3 + p2p4)*
     &                    3D0 + ame**4*8D0 + p1p2**2*8D0 - 
     &                   p1p2*1D0*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     ame**2*1.6D1))) - 
     &             kp2*1D0*
     &              (p1p2**3*
     &                 (p2p4*2D0 + ampi**2*3D0 - p3p4*3D0) + 
     &                p1p2**2*
     &                 (p2p4*(p1p3 + ame**2*2D0)*3D0 + 
     &                   p1p4*(p2p3 - (ame**2 + p2p4)*2D0)*
     &                    3D0 + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*4D0) + 
     &                ame**2*p1p2*
     &                 (p2p4**2*4D0 + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*5D0 - 
     &                   p1p4*2D0*
     &                    (ame**2*2D0 + p1p3*3D0 - p2p3*3D0 + 
     &                     p2p4*5D0) + p1p4**2*6D0 + 
     &                   p2p4*
     &                    (p1p3*4D0 - p2p3*4D0 + ame**2*6D0))
     &                 + ame**4*
     &                 (p1p4**2*2D0 + 
     &                   ame**2*(ampi**2 - p3p4*1D0)*2D0 - 
     &                   p2p4**2*4D0 + 
     &                   p2p4*
     &                    (ame**2*2D0 + p2p3*4D0 - p1p3*7D0)
     &                    + p1p4*
     &                    (-(p2p3*1D0) + ame**2*2D0 - 
     &                     p1p3*2D0 + p2p4*8D0)) + 
     &                kp4*
     &                 (ame**4*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     ame**2*2D0) - 
     &                   ame**2*p1p2*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     ame**2*2D0)*5D0 - p1p2**3*6D0 + 
     &                   ame**2*p1p2**2*3.D1)))) + 
     &       kp1*(kp2**4*
     &           (kp4*(p1p3*(p1p2 - ame**2*1D0) + 
     &                ame**2*(p3p4 - ampi**2*1D0 - p2p3*2D0))
     &              + kp3*
     &              (p1p4*(p1p2 - ame**2*1D0) + 
     &                ame**2*
     &                 (p3p4 - ampi**2*1D0 - p2p4*2D0) + 
     &                ame**2*kp4*4D0)) + 
     &          ame**2*(ame**2 + p1p2)**2*(ame**2 - p1p2*1D0)*
     &           (-(kp3**2*(p1p4 + p2p4)*1D0) - 
     &             kp4*1D0*
     &              (kp4*(p1p3 + p2p3) + 
     &                p2p3*
     &                 (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*2D0
     &                ) + 
     &             kp3*(p2p4*
     &                 (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*2D0
     &                  + kp4*
     &                 (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                   ame**2*4D0 - p1p2*4D0))) - 
     &          kp2*(ame**2 + p1p2)*2D0*
     &           (-(ame**2*kp3**2*(p1p4 + p2p4)*1D0*
     &                (ame**2 - p1p2*3D0)) + 
     &             kp4*(p1p2**3*(ampi**2 - p3p4*1D0) + 
     &                p1p2**2*
     &                 (p2p3*(p1p4 - ame**2*1D0) + 
     &                   p1p3*(p2p4 - ame**2*1D0 - p2p3*2D0))
     &                 + ame**2*p1p2*
     &                 (ame**2*(p3p4 - ampi**2*1D0) + 
     &                   kp4*(p1p3 + p2p3)*3D0 + 
     &                   (p1p3 - p2p3*1D0)*
     &                    (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &                    3D0) + 
     &                ame**4*
     &                 (-(p1p3**2*1D0) + 
     &                   p2p3*
     &                    (ame**2 + p2p4 - kp4*1D0 - 
     &                     p2p3*1D0 - p1p4*2D0) + 
     &                   p1p3*
     &                    (ame**2 + p1p4 - kp4*1D0 - 
     &                     p2p4*2D0 + p2p3*4D0))) + 
     &             kp3*(p1p2**3*(ampi**2 - p3p4*1D0) + 
     &                p1p2**2*
     &                 (p1p4*(p2p3 - ame**2*1D0) + 
     &                   p2p4*(p1p3 - ame**2*1D0 - p1p4*2D0))
     &                 + ame**2*p1p2*
     &                 (ame**2*(p3p4 - ampi**2*1D0) + 
     &                   (p1p4 - p2p4*1D0)*
     &                    (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*
     &                    3D0) + 
     &                ame**4*
     &                 (-(p2p4**2*1D0) + 
     &                   p1p4*
     &                    (ame**2 + p1p3 - p1p4*1D0 - 
     &                     p2p3*2D0) + 
     &                   p2p4*
     &                    (ame**2 + p2p3 - p1p3*2D0 + 
     &                     p1p4*4D0)) + 
     &                kp4*
     &                 (-(p1p2**3*2D0) + 
     &                   ame**4*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     ame**2*6D0) + 
     &                   ame**2*p1p2**2*1.D1 - 
     &                   ame**2*p1p2*1D0*
     &                    (p1p3*3D0 + p1p4*3D0 + p2p3*3D0 + 
     &                     p2p4*3D0 + ame**2*1.4D1)))) - 
     &          kp2**3*1D0*
     &           (ame**2*kp3**2*(p1p4 + p2p4) + 
     &             kp4*(p1p2**2*
     &                 (ampi**2 - p3p4*1D0 + p1p3*3D0) + 
     &                ame**2*
     &                 (p2p3**2*2D0 + 
     &                   ame**2*(ampi**2 - p3p4*1D0)*2D0 + 
     &                   p2p3*
     &                    (kp4 + p1p4 - (ame**2 + p2p4)*2D0)
     &                    + p1p3*(kp4 - p2p4*1D0 - ame**2*3D0)
     &                   ) + 
     &                p1p2*
     &                 (p1p3*(p2p4 - p2p3*2D0 + ame**2*4D0) + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*5D0 + 
     &                   p2p3*(p1p4 - ame**2*6D0))) + 
     &             kp3*(p1p2**2*
     &                 (ampi**2 - p3p4*1D0 + p1p4*3D0) + 
     &                ame**2*
     &                 (p2p4**2*2D0 + 
     &                   ame**2*(ampi**2 - p3p4*1D0)*2D0 + 
     &                   p2p4*(p1p3 - (ame**2 + p2p3)*2D0) - 
     &                   p1p4*1D0*(p2p3 + ame**2*3D0)) + 
     &                p1p2*
     &                 (p1p4*(p2p3 - p2p4*2D0 + ame**2*4D0) + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*5D0 + 
     &                   p2p4*(p1p3 - ame**2*6D0)) - 
     &                kp4*1D0*
     &                 (p1p2**2*2D0 + 
     &                   ame**2*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     ame**2*1.D1) - ame**2*p1p2*2.D1)))
     &           + kp2**2*
     &           (-(ame**2*kp3**2*(p1p4 + p2p4)*1D0*
     &                (ame**2 - p1p2*5D0)) + 
     &             kp4*(p1p2**3*
     &                 (p1p3*2D0 + ampi**2*3D0 - p3p4*3D0) + 
     &                p1p2**2*
     &                 (p2p3*(p1p4 - ame**2*2D0)*3D0 + 
     &                   p1p3*(p2p4 + ame**2*2D0 - p2p3*2D0)*
     &                    3D0 + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*4D0) + 
     &                ame**4*
     &                 (ame**2*ampi**2*2D0 - 
     &                   ame**2*p3p4*2D0 + 
     &                   p2p3*
     &                    (-(kp4*1D0) - p1p4*1D0 + 
     &                     (ame**2 + p2p3 - p2p4*1D0)*2D0) - 
     &                   p1p3**2*4D0 + 
     &                   p1p3*
     &                    (-(kp4*1D0) + ame**2*2D0 + 
     &                     p1p4*4D0 - p2p4*7D0 + p2p3*8D0)) + 
     &                ame**2*p1p2*
     &                 (p1p3**2*4D0 + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*5D0 + 
     &                   p2p3**2*6D0 + 
     &                   p2p3*
     &                    (-(ame**2*4D0) + kp4*5D0 + 
     &                     p1p4*6D0 - p2p4*6D0) + 
     &                   p1p3*
     &                    (-(p1p4*4D0) + p2p4*4D0 + kp4*5D0 + 
     &                     ame**2*6D0 - p2p3*1.D1))) + 
     &             kp3*(p1p2**3*
     &                 (p1p4*2D0 + ampi**2*3D0 - p3p4*3D0) + 
     &                p1p2**2*
     &                 (p2p4*(p1p3 - ame**2*2D0)*3D0 + 
     &                   p1p4*(p2p3 + ame**2*2D0 - p2p4*2D0)*
     &                    3D0 + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*4D0) + 
     &                ame**2*p1p2*
     &                 (p1p4**2*4D0 + 
     &                   ame**2*(p3p4 - ampi**2*1D0)*5D0 + 
     &                   p1p4*2D0*
     &                    (-(p1p3*2D0) + p2p3*2D0 + 
     &                     ame**2*3D0 - p2p4*5D0) + 
     &                   p2p4**2*6D0 + 
     &                   p2p4*
     &                    (-(ame**2*4D0) + p1p3*6D0 - 
     &                     p2p3*6D0)) + 
     &                ame**4*
     &                 (ame**2*ampi**2*2D0 - 
     &                   ame**2*p3p4*2D0 + 
     &                   p2p4*
     &                    (-(p1p3*1D0) + 
     &                     (ame**2 + p2p4 - p2p3*1D0)*2D0) - 
     &                   p1p4**2*4D0 + 
     &                   p1p4*
     &                    (ame**2*2D0 + p1p3*4D0 - p2p3*7D0 + 
     &                     p2p4*8D0)) + 
     &                kp4*
     &                 (ame**4*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     ame**2*2D0) - 
     &                   ame**2*p1p2*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     ame**2*2D0)*5D0 - p1p2**3*6D0 + 
     &                   ame**2*p1p2**2*3.D1))))))/
     &   (kp1*kp2*kp3*kp4*(ame**2 + p1p2)*
     &     (-(ame**4*1D0) + (kp2 - p1p2*1D0)**2 + 
     &       ame**2*kp1*2D0)*
     &     (-(ame**4*1D0) + (kp1 - p1p2*1D0)**2 + 
     &       ame**2*kp2*2D0)) - 
     &  (C0aex25x34m1em2*ep2*5.D-1*
     &     (ame**2*kp2**3*(kp4*p1p3 + kp3*p1p4)*4D0 + 
     &       ame**2*kp1**3*4D0*
     &        (kp4*(p3p4 - ampi**2*1D0 + p2p3*2D0) + 
     &          kp3*(p3p4 - ampi**2*1D0 + p2p4*2D0 + kp4*4D0))
     &         + ame**2*(ame**2 + p1p2)*
     &        (m12 + m22 - p1p2*4D0)*
     &        (kp3**2*(p1p4 + p2p4) + 
     &          kp4*(kp4*(p1p3 + p2p3) + 
     &             p2p3*(p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*
     &              2D0) + 
     &          kp3*(p2p4*(p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &              2D0 - 
     &             kp4*1D0*
     &              (p1p3 + p1p4 + p2p3 + p2p4 + ame**2*4D0 - 
     &                p1p2*4D0))) - 
     &       kp2**2*2D0*(kp4*
     &           (p1p2*(p1p3*(-(m12*1D0) + ame**2*2D0) + 
     &                ame**2*(ampi**2 - p3p4*1D0 - p2p3*2D0))
     &              + ame**2*
     &              ((ame**2 + m12 + p1p4)*p2p3*2D0 + 
     &                (ampi**2 - p3p4*1D0)*
     &                 (m12 + ame**2*3D0) + 
     &                p1p3*
     &                 (-(m22*1D0) + ame**2*2D0 + p2p4*2D0 - 
     &                   p2p3*4D0))) + 
     &          kp3*(p1p2*
     &              (p1p4*(-(m12*1D0) + ame**2*2D0) + 
     &                ame**2*(ampi**2 - p3p4*1D0 - p2p4*2D0))
     &              - ame**2*kp4*4D0*(m12 + ame**2*4D0) + 
     &             ame**2*
     &              ((ame**2 + m12 + p1p3)*p2p4*2D0 + 
     &                (ampi**2 - p3p4*1D0)*
     &                 (m12 + ame**2*3D0) + 
     &                p1p4*
     &                 (-(m22*1D0) + ame**2*2D0 + p2p3*2D0 - 
     &                   p2p4*4D0)))) + 
     &       kp2*(-(ame**2*kp3**2*(p1p4 + p2p4)*
     &             (ame**2 + m12 - p1p2*1D0)*2D0) + 
     &          kp4*(ame**2*2D0*
     &              (-(kp4*(ame**2 + m12)*(p1p3 + p2p3)*
     &                   1D0) + 
     &                p2p4*
     &                 ((ame**2 + m12)*p2p3*2D0 + 
     &                   p1p3*(-(m22*1D0) + ame**2*2D0)) + 
     &                p2p3*
     &                 (ame**2*(m12 + m22) - 
     &                   (ame**2 + m12)*p2p3*2D0 + 
     &                   p1p3*(m12 + m22 - ame**2*1D0)*2D0 - 
     &                   p1p4*1D0*(m22 + m12*2D0))) + 
     &             ame**4*p3p4*(m22 - m12*1D0 - ame**2*4D0) + 
     &             ame**4*ampi**2*
     &              (m12 - m22*1D0 + ame**2*4D0) - 
     &             p1p2**2*2D0*
     &              (m12*(ampi**2 - p3p4*1D0) + 
     &                ame**2*p2p3*4D0) + 
     &             p1p2*(ame**2*kp4*(p1p3 + p2p3)*2D0 - 
     &                p2p4*2D0*
     &                 (ame**2*p2p3*2D0 + 
     &                   p1p3*(m12 - ame**2*2D0)) + 
     &                ame**2*p3p4*(m12 + m22 - ame**2*4D0) + 
     &                ame**2*ampi**2*
     &                 (-(m12*1D0) - m22*1D0 + ame**2*4D0) + 
     &                p2p3*2D0*
     &                 (p1p3*2D0*(m12 - ame**2*3D0) + 
     &                   ame**2*
     &                    (m12 + m22 + p2p3*2D0 - ame**2*4D0)
     &                    + p1p4*(-(m12*1D0) + ame**2*4D0))))
     &           + kp3*(-(p1p2**2*2D0*
     &                (m12*(ampi**2 - p3p4*1D0) + 
     &                  ame**2*p2p4*4D0)) + 
     &             p1p2*(p1p4*2D0*
     &                 (p2p3*(-(m12*1D0) + ame**2*2D0) + 
     &                   p2p4*2D0*(m12 - ame**2*3D0)) + 
     &                ame**2*p2p4**2*4D0 + 
     &                ame**2*(ampi**2 - p3p4*1D0)*
     &                 (-(m12*1D0) - m22*1D0 + ame**2*4D0) + 
     &                p2p4*2D0*
     &                 (ame**2*
     &                    (m12 + m22 - p2p3*2D0 - ame**2*4D0)
     &                    + p1p3*(-(m12*1D0) + ame**2*4D0)))
     &              + ame**2*
     &              (p2p4*2D0*
     &                 (ame**2*(m12 + m22) + 
     &                   (ame**2 + m12)*p2p3*2D0 - 
     &                   p1p3*1D0*(m22 + m12*2D0)) - 
     &                (ame**2 + m12)*p2p4**2*4D0 + 
     &                ame**2*(ampi**2 - p3p4*1D0)*
     &                 (m12 - m22*1D0 + ame**2*4D0) + 
     &                p1p4*
     &                 (p2p4*(m12 + m22 - ame**2*1D0)*4D0 + 
     &                   p2p3*(-(m22*2D0) + ame**2*4D0))) + 
     &             kp4*(ame**2*(ame**2 + m12)*
     &                 (p1p3 + p1p4 + p2p3 + p2p4)*2D0 + 
     &                ame**4*4D0*
     &                 (m12 + m22*2D0 - ame**2*4D0) + 
     &                p1p2**2*4D0*(m12 + ame**2*4D0) - 
     &                ame**2*p1p2*2D0*
     &                 (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                   ame**2*1.6D1)))) + 
     &       kp1**2*(-(ame**2*kp3**2*(p1p4 + p2p4)*4D0) + 
     &          kp4*(-(kp2*2D0*
     &                (m12*p2p3 + 
     &                  ame**2*(ampi**2 - p3p4*1D0)*3D0)) + 
     &             ame**2*
     &              (-(kp4*(p1p3 + p2p3)*4D0) + 
     &                (ampi**2 - p3p4*1D0 - p2p3*2D0)*
     &                 (-(m12*1D0) - m22*1D0 + ame**2*4D0) + 
     &                p2p3*
     &                 (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*8D0
     &                  + p1p2*
     &                 (ampi**2 - p3p4*1D0 - p2p3*2D0)*8D0))
     &           + kp3*(ame**2*
     &              ((ampi**2 - p3p4*1D0 - p2p4*2D0)*
     &                 (-(m12*1D0) - m22*1D0 + ame**2*4D0) + 
     &                p2p4*
     &                 (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*8D0
     &                  + p1p2*
     &                 (ampi**2 - p3p4*1D0 - p2p4*2D0)*8D0 + 
     &                kp4*4D0*
     &                 (m12 + m22 + p1p3 + p1p4 + p2p3 + 
     &                   p2p4 - p1p2*1.2D1)) + 
     &             kp2*(-(m12*p2p4*2D0) + 
     &                ame**2*(p3p4 - ampi**2*1D0)*6D0 + 
     &                ame**2*kp4*2.4D1))) + 
     &       kp1*(ame**2*kp3**2*(p1p4 + p2p4)*
     &           (-(m12*1D0) - m22*1D0 - kp2*2D0 + 
     &             ame**2*4D0 + p1p2*8D0) + 
     &          kp4*(-(kp2**2*2D0*
     &                (ame**2*(ampi**2 - p3p4*1D0) + 
     &                  p1p3*(m12 - ame**2*2D0))) + 
     &             kp2*(-(ame**2*kp4*(p1p3 + p2p3)*2D0) + 
     &                ame**2*p3p4*(m12 + m22 + ame**2*2D0) - 
     &                ame**2*ampi**2*1D0*
     &                 (m12 + m22 + ame**2*2D0) + 
     &                p2p4*2D0*
     &                 (ame**2*p2p3*2D0 + 
     &                   p1p3*(m12 - ame**2*2D0)) + 
     &                p1p2*2D0*
     &                 ((ampi**2 - p3p4*1D0)*
     &                    (m12 + ame**2*3D0) + 
     &                   p2p3*(m12 + ame**2*4D0)) + 
     &                p2p3*2D0*
     &                 (p1p4*(m12 - ame**2*4D0) + 
     &                   ame**2*
     &                    (m12 - p2p3*2D0 + ame**2*4D0) + 
     &                   p1p3*(-(m12*2D0) + ame**2*6D0))) + 
     &             ame**2*
     &              (ame**2*ampi**2*(m12 + m22) - 
     &                ame**2*(m12 + m22)*p3p4*1D0 + 
     &                p1p2**2*(p3p4 - ampi**2*1D0 + p2p3*2D0)*
     &                 4D0 + 
     &                kp4*(p1p3 + p2p3)*
     &                 (-(m12*1D0) - m22*1D0 + ame**2*4D0) - 
     &                p2p3*2D0*
     &                 (ame**2*(m12 + m22) + 
     &                   (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &                    (-(m12*1D0) - m22*1D0 + ame**2*4D0))
     &                  + p1p2*
     &                 (-(1D0*(ampi**2 - p3p4*1D0 - p2p3*2D0)*
     &                     (-(m12*1D0) - m22*1D0 + ame**2*4D0)
     &                     ) + kp4*(p1p3 + p2p3)*8D0 + 
     &                   p2p3*
     &                    (p1p4 + p2p3 - p1p3*1D0 - p2p4*1D0)*
     &                    1.6D1))) + 
     &          kp3*(kp2**2*
     &              (ame**2*(p3p4 - ampi**2*1D0)*2D0 + 
     &                p1p4*(-(m12*2D0) + ame**2*4D0) + 
     &                ame**2*kp4*8D0) + 
     &             kp2*(p1p4*p2p3*2D0*(m12 - ame**2*2D0) + 
     &                ame**2*p3p4*(m12 + m22 + ame**2*2D0) - 
     &                ame**2*ampi**2*1D0*
     &                 (m12 + m22 + ame**2*2D0) - 
     &                ame**2*p2p4**2*4D0 + 
     &                p1p2*2D0*
     &                 ((ampi**2 - p3p4*1D0)*
     &                    (m12 + ame**2*3D0) + 
     &                   p2p4*(m12 + ame**2*4D0)) + 
     &                p2p4*2D0*
     &                 (p1p3*(m12 - ame**2*4D0) + 
     &                   ame**2*
     &                    (m12 + p2p3*2D0 + ame**2*4D0) + 
     &                   p1p4*(-(m12*2D0) + ame**2*6D0)) + 
     &                kp4*2D0*
     &                 (-(p1p2*2D0*(m12 + ame**2*1.D1)) + 
     &                   ame**2*
     &                    (p1p3 + p1p4 + p2p3 + p2p4 + 
     &                     m22*2D0 + ame**2*1.2D1))) + 
     &             ame**2*
     &              (ame**2*ampi**2*(m12 + m22) - 
     &                ame**2*(m12 + m22)*p3p4*1D0 + 
     &                p1p2**2*(p3p4 - ampi**2*1D0 + p2p4*2D0)*
     &                 4D0 + 
     &                p2p4*2D0*
     &                 (-(ame**2*(m12 + m22)*1D0) + 
     &                   (p1p3 + p2p4 - p1p4*1D0 - p2p3*1D0)*
     &                    (-(m12*1D0) - m22*1D0 + ame**2*4D0))
     &                  + p1p2*
     &                 (-(1D0*(ampi**2 - p3p4*1D0)*
     &                     (-(m12*1D0) - m22*1D0 + ame**2*4D0)
     &                     ) + 
     &                   p2p4*2D0*
     &                    (-(m12*1D0) - m22*1D0 + 
     &                     ame**2*4D0 + p1p3*8D0 - p1p4*8D0 - 
     &                     p2p3*8D0) + p2p4**2*1.6D1) + 
     &                kp4*
     &                 (-((p1p3 + p1p4 + p2p3 + p2p4)*1D0*
     &                     (-(m12*1D0) - m22*1D0 + ame**2*4D0)
     &                     ) - 
     &                   p1p2*
     &                    (m12 + m22 + p1p3 + p1p4 + p2p3 + 
     &                     p2p4)*8D0 - ame**4*1.6D1 + 
     &                   p1p2**2*4.8D1))))))/
     &   (kp2*kp3*kp4*(ame**2 + p1p2)*
     &     (-(ame**4*1D0) + (kp1 - p1p2*1D0)**2 + 
     &       ame**2*kp2*2D0))
                  
