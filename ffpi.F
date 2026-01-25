      !!==========================!!
      !!                          !!
      !!   FORM FACTOR ROUTINES   !!
      !!                          !!
      !!==========================!!

      !!==============================================================!!

      !! Charged pion form factor
      !! Input: squared momentum (q2)
      !! Output: complex form factor (ff), squared form factor (ff2)
      !! Implemented parametrisations: strong2020, cmd2, cmd3, bwsum2,
      !! bwsum3, snd, babar, besiii, kloe2, phokhara, bern
      subroutine ffpi(q2,ff,ff2)

      implicit none

      double precision q2,ff2,ff2pr,q2pr
      double complex ff,ffpr

      double precision toll
      parameter (toll = 1d-10)

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi

      common/previousff/q2pr,ffpr,ff2pr
      data q2pr,ffpr,ff2pr/-1.d10,(1.d0,0.d0),1.d0/

      ff  = (1.d0,0.d0)
      ff2 = 1.d0

      !! No pion form factor
      if (iffpi.eq.0) return

      if (abs((q2-q2pr)/q2).gt.toll) then

      !! Select the pion form factor parametrisation
         if (what_ffpi.eq.'strong2020') then
            call ffpi_strong2020(q2,ff,ff2)
         elseif (what_ffpi.eq.'cmd2') then
            call ffpi_cmd2(q2,ff,ff2)
         elseif (what_ffpi.eq.'cmd3') then
            call ffpi_cmd3(q2,ff,ff2)
         elseif (what_ffpi.eq.'bwsum2') then
            call ffpi_bwsum2(q2,ff,ff2)
         elseif (what_ffpi.eq.'bwsum3') then
            call ffpi_bwsum3(q2,ff,ff2)
         elseif (what_ffpi.eq.'snd') then
            call ffpi_snd(q2,ff,ff2)
         elseif (what_ffpi.eq.'babar') then
            call ffpi_babar(q2,ff,ff2)
         elseif (what_ffpi.eq.'besiii') then
            call ffpi_besiii(q2,ff,ff2)
         elseif (what_ffpi.eq.'kloe2') then
            call ffpi_kloe2(q2,ff,ff2)
         elseif (what_ffpi.eq.'phokhara') then
            call ffpi_phokhara(q2,ff,ff2)
         elseif (what_ffpi.eq.'bern') then
            call ffpi_bern(q2,ff,ff2)
         else
            print*, 'Error: Invalid choice of pion form factor'
         endif

         q2pr  = q2
         ffpr  = ff
         ff2pr = ff2

      else

         ff  = ffpr
         ff2 = ff2pr

      endif

      return
      end

      !!==============================================================!!

      !! Charged pion form factor as in PHOKHARA 10.0
      !! Data from CMD-2, BABAR, KLOE-2, SND, DM2, CLEO-c
      !! Reference: Phys.Rev.D 81 (2010) 094014
      !! Input: squared momentum (q2)
      !! Output: complex form factor (ff), squared form factor (ff2)
      !! Energy range: 4*m_pi^2 < s < 16 GeV^2
      subroutine ffpi_phokhara(q2,ff,ff2)

      implicit none

      double precision q2,ff2
      double complex ff

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      double precision ampi
      common/hparam/ampi

      double precision beta,iloc
      common/par_dqcd/beta
      integer idqcd,ndqcd

      double precision m_n(0:5),g_n(0:5)
      double precision mod_n(0:5),arg_n(0:5),c_n(0:5)
      double complex f_n(0:5)

      double complex bw_const,bw_gs,ff_dqcd

      double precision pi
      integer k,ifirst

      double precision s_max
      common/smax/s_max

      common/ffpiphokharacmn/pi,m_n,g_n,mod_n,arg_n,f_n,c_n,ifirst
      data ifirst /0/

      !! Maximum safe q^2 (GeV^2)
      s_max = 16.d0
      call ffpi_warning(q2,s_max)

      !! Flag to activate the contribution from Dual QCD for n >= 6
      !! It can be computationally expensive
      idqcd = 1

      !! Number of resonances to be computed in Dual QCD
      !! It must be ndqcd >= 6 (in PHOKHARA ndqcd = 2000)
      ndqcd = 2000

      if (ifirst.eq.0) then

         pi = 4.d0*atan(1.d0)

         !! Charged pion mass (GeV)
         m_pi = ampi

         !! Omega mass and width (GeV)
         m_om = 0.78241d0
         g_om = 0.00833d0

         !! Dual QCD input parameter from fit
         beta = 2.14845339006279d0

         !! Masses of vector mesons (GeV)
         m_n(0) = 0.7733690734480345d0
         m_n(1) = 1.489505729274754d0
         m_n(2) = 1.86961155147573d0
         m_n(3) = 2.12d0
         m_n(4) = 2.320153904230950d0
         m_n(5) = 2.565026652070740d0

         !! Widths of vector mesons(GeV)
         g_n(0) = 0.1470530254686560d0
         g_n(1) = 0.4290884148991568d0
         g_n(2) = 0.3575024923053209d0
         g_n(3) = 0.3d0
         g_n(4) = 0.4414835430304250d0
         g_n(5) = 0.4880784211162140d0

         !! Rho mass and width (GeV)
         m_rho = m_n(0)
         g_rho = g_n(0)

         !! Resonance amplitudes

         mod_n(1) = 0.001867242111672990d0
         mod_n(2) = 0.5894475353098135d0
         mod_n(3) = 4.767753409025326d-2
         mod_n(4) = 0.3962340510662401d0
         mod_n(5) = 0.4302294935131883d0

         arg_n(1) = 0.1055961852734719d0
         arg_n(2) = -2.199783053184180d0
         arg_n(3) = -2.040841165327912d0
         arg_n(4) = -2.944688084027720d0
         arg_n(5) = 1.191141958695573d0

         f_n(1) = mod_n(1)*dcmplx(cos(arg_n(1)),sin(arg_n(1)))
         f_n(2) = mod_n(2)*dcmplx(cos(arg_n(2)),sin(arg_n(2)))
         f_n(3) = mod_n(3)*dcmplx(cos(arg_n(3)),sin(arg_n(3)))
         f_n(4) = mod_n(4)*dcmplx(cos(arg_n(4)),sin(arg_n(4)))
         f_n(5) = mod_n(5)*dcmplx(cos(arg_n(5)),sin(arg_n(5)))

         !! Resonance amplitudes from Dual QCD

         iloc = 1.d0
         c_n(0) = 1.087633403691967d0

         do k=1,5
            iloc = iloc * (1.d0 - (beta - 1.d0)/k)
            c_n(k) = (-1.d0)**k * gamma(beta - 0.5d0) * 2.d0/sqrt(pi)
     &           / (1.d0 + 2.d0*k) * iloc/pi * gamma(2.d0 - beta)
     &           * sin(pi*(beta - 1.d0 - k))
         enddo
         ifirst = 1
      endif

      !! Pion form factor
      ff = c_n(0) * bw_gs(m_n(0),g_n(0),q2) / (1.d0 + f_n(1))
     &   * (1.d0 + f_n(1) * bw_const(m_om,g_om,q2))
     &   + (c_n(1) + c_n(2) + c_n(3) + c_n(4) + c_n(5))
     &   / (1.d0 + f_n(2) + f_n(3) + f_n(4) + f_n(5))
     &   * (bw_gs(m_n(1),g_n(1),q2)
     &   + f_n(2) * bw_gs(m_n(2),g_n(2),q2)
     &   + f_n(3) * bw_gs(m_n(3),g_n(3),q2)
     &   + f_n(4) * bw_gs(m_n(4),g_n(4),q2)
     &   + f_n(5) * bw_gs(m_n(5),g_n(5),q2))

      !! Add Dual QCD contribution if activated
      if (idqcd.eq.1.and.ndqcd.gt.5) then
         ff = ff + ff_dqcd(q2,ndqcd)
      endif

      !! Squared pion form factor
      ff2 = real(ff)**2 + aimag(ff)**2

      return
      end

      !!==============================================================!!

      !! Charged pion form factor from CMD-2
      !! Reference: Phys.Lett.B 648 (2007) 28-38,
      !!            JETP Lett. 82 (2005) 743-747
      !! Input: squared momentum (q2)
      !! Output: complex form factor (ff), squared form factor (ff2)
      !! Energy range: 4*m_pi^2 < s < 2 GeV^2
      subroutine ffpi_cmd2(q2,ff,ff2)

      implicit none

      double precision q2,ff2
      double complex ff

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      double precision ampi
      common/hparam/ampi

      double precision m_n(0:1),g_n(0:1)
      double precision mod_n(0:1),arg_n(0:1)
      double complex f_n(0:1)

      double precision mod_om,arg_om,beta
      double complex del_om

      double complex bw_const,bw_gs

      double precision m_phi,g_phi
      double precision mod_phi,arg_phi
      double complex del_phi

      double precision s_max
      common/smax/s_max

      !! Maximum safe q^2 (GeV^2)
      s_max = 2.d0
      call ffpi_warning(q2,s_max)

      !! Charged pion mass (GeV)
      m_pi = ampi

      !! Omega mass and width (GeV)
      m_om = 0.782659615158d0
      g_om = 0.008549875257d0

      !! Masses of vector mesons (GeV)
      m_n(0) = 0.7754924021d0
      m_n(1) = 1.3697607590d0

      !! Widths of vector mesons (GeV)
      g_n(0) = 0.1456898531d0
      g_n(1) = 0.3852126547d0

      !! Resonance amplitudes

      mod_n(1) = 0.08767160d0
      arg_n(1) = 3.15866288d0

      mod_om = 0.00157628957d0
      arg_om = 0.17877329107d0

      f_n(1) = mod_n(1)*dcmplx(cos(arg_n(1)),sin(arg_n(1)))
      del_om = mod_om*dcmplx(cos(arg_om),sin(arg_om))

      !! Pion form factor
      ff = (bw_gs(m_n(0),g_n(0),q2)
     &   * (1 + del_om*q2/m_om**2 * bw_const(m_om,g_om,q2))
     &   + f_n(1) * bw_gs(m_n(1),g_n(1),q2)) / (1 + f_n(1))

      !! Squared pion form factor
      ff2 = real(ff)**2 + aimag(ff)**2

      return
      end

      !!==============================================================!!

      !! Charged pion form factor from CMD-3 (+ CMD-2 and DM2)
      !! Reference: ArXiv:2302.08834, arXiv:2309.12910
      !! Input: squared momentum (q2)
      !! Output: complex form factor (ff), squared form factor (ff2)
      !! Energy range: 4*m_pi^2 < s < 4 GeV^2
      subroutine ffpi_cmd3(q2,ff,ff2)

      implicit none

      double precision q2,ff2
      double complex ff

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      double precision ampi
      common/hparam/ampi

      double precision m_n(0:2),g_n(0:2)
      double precision mod_n(0:2),arg_n(0:2)
      double complex f_n(0:2)

      double precision m_phi,g_phi
      double precision mod_om,arg_om
      double precision mod_phi,arg_phi
      double complex del_om,del_phi

      double complex bw_const,bw_gs

      double precision s_max
      common/smax/s_max

      !! Maximum safe q^2 (GeV^2)
      s_max = 4.d0
      call ffpi_warning(q2,s_max)

      !! Charged pion mass (GeV)
      m_pi = ampi

      !! Omega mass and width (GeV)
      m_om = 0.782239642d0
      g_om = 0.008179181218d0

      !! Phi mass and width (GeV)
      m_phi = 1.0195515864d0
      g_phi = 0.005361127785d0

      !! Masses of vector mesons (GeV)
      m_n(0) = 0.773947746d0
      m_n(1) = 1.4566860687d0
      m_n(2) = 1.8707382947d0

      !! Widths of vector mesons (GeV)
      g_n(0) = 0.14809686121d0
      g_n(1) = 0.5240547792d0
      g_n(2) = 0.17048893700d0

      !! Resonance amplitudes

      mod_n(1) = 0.09785638624d0
      mod_n(2) = 0.03744531577d0
      mod_om = 0.001643118856d0
      mod_phi = 0.00059360578740d0

      arg_n(1) = 3.5412996691092d0
      arg_n(2) = 2.2770437268d0
      arg_om = 0.05715639453d0
      arg_phi = 2.836459685d0

      f_n(0) = (1.d0,0.d0)
      f_n(1) = mod_n(1)*dcmplx(cos(arg_n(1)),sin(arg_n(1)))
      f_n(2) = mod_n(2)*dcmplx(cos(arg_n(2)),sin(arg_n(2)))
      del_om = mod_om*dcmplx(cos(arg_om),sin(arg_om))
      del_phi = mod_phi*dcmplx(cos(arg_phi),sin(arg_phi))

      !! Pion form factor
      ff = (bw_gs(m_n(0),g_n(0),q2)
     &   * (f_n(0) + del_om*q2/m_om**2 * bw_const(m_om,g_om,q2)
     &   + del_phi*q2/m_phi**2 * bw_const(m_phi,g_phi,q2))
     &   + f_n(1) * bw_gs(m_n(1),g_n(1),q2)
     &   + f_n(2) * bw_gs(m_n(2),g_n(2),q2))
     &   / (f_n(0) + f_n(1) + f_n(2))

      !! Squared pion form factor
      ff2 = real(ff)**2 + aimag(ff)**2

      return
      end

      !!==============================================================!!

      !! Fit of the CMD-2 form factor as a sum of simple BW functions
      !! To be used in the GVMD approach
      !! Input: squared momentum (q2)
      !! Output: complex form factor (ff), squared form factor (ff2)
      !! Energy range: 4*m_pi^2 < s < 2 GeV^2
      subroutine ffpi_bwsum2(q2,ff,ff2)

      implicit none

      double precision q2,ff2
      double complex ff

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      double precision ampi
      common/hparam/ampi

      integer i,nres
      double precision m_n(1:9),g_n(1:9)
      double precision mod_n(1:9),arg_n(1:9)
      double complex a_n(1:9),a_sum
      common/par_gvmd/m_n,g_n,a_n

      double complex bw_gvmd

      double precision s_max
      common/smax/s_max

      double complex im
      im = (0.d0,1.d0)

      !! Maximum safe q^2 (GeV^2)
      s_max = 2.d0
      call ffpi_warning(q2,s_max)

      !! Charged pion mass (GeV)
      m_pi = ampi

      !! Number of resonances
      nres = 3

      !! Masses of vector mesons (GeV)
      m_n(1) = 0.7580789340642422d0
      m_n(2) = 0.782802909129871d0
      m_n(3) = 1.2538111054640864d0

      !! Widths of vector mesons (GeV)
      g_n(1) = 0.1368056759510071d0
      g_n(2) = 0.008003969603568873d0
      g_n(3) = 0.5308649685845045d0

      !! Resonance amplitudes
      mod_n(2) = 0.007944883506779365d0
      mod_n(3) = 0.14410469499004908d0

      arg_n(2) = 2.0142117794789183d0
      arg_n(3) = 3.0207621669132627d0

      a_n(1) = (1.d0,0.d0)
      a_sum = a_n(1)
      do i=2,nres
         a_n(i) = mod_n(i)*exp(im*arg_n(i))
         a_sum = a_sum + a_n(i)
      enddo

      !! Form factor
      ff = (0.d0, 0.d0)
      do i=1,nres
         a_n(i) = a_n(i)/a_sum
         ff = ff + a_n(i)*bw_gvmd(m_n(i),g_n(i),q2)
      enddo

      !! Squared pion form factor
      ff2 = real(ff)**2 + aimag(ff)**2

      return
      end

      !!==============================================================!!

      !! Fit of the CMD-3 form factor as a sum of simple BW functions
      !! To be used in the GVMD approach
      !! Input: squared momentum (q2)
      !! Output: complex form factor (ff), squared form factor (ff2)
      !! Energy range: 4*m_pi^2 < s < 2 GeV^2
      subroutine ffpi_bwsum3(q2,ff,ff2)

      implicit none

      double precision q2,ff2
      double complex ff

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      double precision ampi
      common/hparam/ampi

      integer i,nres
      double precision m_n(1:9),g_n(1:9)
      double precision mod_n(1:9),arg_n(1:9)
      double complex a_n(1:9),a_sum
      common/par_gvmd/m_n,g_n,a_n

      double complex bw_gvmd

      double precision s_max
      common/smax/s_max

      double complex bw_const

      double complex im
      im = (0.d0,1.d0)

      !! Maximum safe q^2 (GeV^2)
      s_max = 2.d0
      call ffpi_warning(q2,s_max)

      !! Charged pion mass (GeV)
      m_pi = ampi

      !! Number of resonances
      nres = 4

      !! Masses of vector mesons (GeV)
      m_n(1) = 0.7562711962277819d0
      m_n(2) = 0.7820318217907292d0
      m_n(3) = 1.243494077338154d0
      m_n(4) = 1.0182278004558887d0

      !! Widths of vector mesons (GeV)
      g_n(1) = 0.13684653236424352d0
      g_n(2) = 0.00798814295403398d0
      g_n(3) = 0.5311217444961865d0
      g_n(4) = 0.00526898216509165d0

      !! Resonance amplitudes

      mod_n(2) = 0.00797813570045449d0
      mod_n(3) = 0.143998137038203d0
      mod_n(4) = 0.000708838282524463d0

      arg_n(2) = 2.0140526209131369d0
      arg_n(3) = 3.01886279574077d0
      arg_n(4) = 5.01886279574077d0

      a_n(1) = (1.d0,0.d0)
      a_sum = a_n(1)
      do i=2,nres
         a_n(i) = mod_n(i)*exp(im*arg_n(i))
         a_sum = a_sum + a_n(i)
      enddo

      !! Form factor
      ff = (0.d0, 0.d0)
      do i=1,nres
         a_n(i) = a_n(i)/a_sum
         ff = ff + a_n(i)*bw_gvmd(m_n(i),g_n(i),q2)
      enddo

      !! Squared pion form factor
      ff2 = real(ff)**2 + aimag(ff)**2

      return
      end

      !!==============================================================!!

      !! Charged pion form factor from SND
      !! Reference: JHEP 01 (2021) 113
      !! Input: squared momentum (q2)
      !! Output: complex form factor (ff), squared form factor (ff2)
      !! Energy range: 4*m_pi^2 < s < 1 GeV^2
      subroutine ffpi_snd(q2,ff,ff2)

      implicit none

      double precision q2,ff2
      double complex ff

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      double precision ampi
      common/hparam/ampi

      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

      double precision m_n(0:1),g_n(0:1),sigma_n(0:1)
      double precision arg_om,sigma_om
      double complex del_om

      double precision q_pi,q_0,q_1,q_om
      double precision g_0_q2,g_1_q2,g_om_q2
      double complex ff_0,ff_1,ff_om

      double precision s_max
      common/smax/s_max

      double complex i

      pi = 4.d0*atan(1.d0)
      i = (0.d0,1.d0)

      !! Maximum safe q^2 (GeV^2)
      s_max = 1.d0
      call ffpi_warning(q2,s_max)

      !! Charged pion mass (GeV)
      m_pi = ampi

      !! Omega mass and width (GeV)
      m_om = 0.78265d0
      g_om = 0.00849d0

      !! Masses of vector mesons (GeV)
      m_n(0) = 0.7753d0
      m_n(1) = 1.465d0

      !! Widths of vector mesons (GeV)
      g_n(0) = 0.1456d0
      g_n(1) = 0.4d0

      !! Resonance amplitudes
      sigma_n(0) = 0.003055375d0
      sigma_n(1) = 6.1636547d-6
      sigma_om = 8.0897968d-5
      arg_om = 1.93207948d0

      del_om = dcmplx(cos(arg_om),sin(arg_om))

      !! Rho contribution

      q_pi = 0.5d0*sqrt(q2 - 4.d0*m_pi**2)
      q_0 = 0.5d0*sqrt(m_n(0)**2 - 4.d0*m_pi**2)
      g_0_q2 = g_n(0)*m_n(0)**2/q2*(q_pi/q_0)**3

      ff_0 = g_n(0)*m_n(0)**3 * sqrt(m_n(0)*sigma_n(0))
     &     / (m_n(0)**2 - q2 - i*sqrt(q2)*g_0_q2) / sqrt(q_0)**3

      !! Rho' contribution

      q_1 = 0.5d0*sqrt(m_n(1)**2 - 4.d0*m_pi**2)
      g_1_q2 = g_n(1)*m_n(1)**2/q2*(q_pi/q_1)**3

      ff_1 = - g_n(1)*m_n(1)**3 * sqrt(m_n(1)*sigma_n(1))
     &     / (m_n(1)**2 - q2 - i*sqrt(q2)*g_1_q2) / sqrt(q_1)**3

      !! Omega contribution

      q_om = 0.5d0*sqrt(m_om**2 - 4.d0*m_pi**2)
      g_om_q2 = g_om*m_om**2/q2*(q_pi/q_om)**3

      ff_om = g_om*m_om**3 * sqrt(m_om*sigma_om) * del_om
     &      / (m_om**2 - q2 - i*sqrt(q2)*g_om_q2) / sqrt(q_om)**3

      !! Pion form factor
      ff = (ff_0 + ff_1 + ff_om) * sqrt(3.d0/8.d0/pi)/alpha

      !! Squared pion form factor
      ff2 = real(ff)**2 + aimag(ff)**2

      return
      end

      !!==============================================================!!

      !! Charged pion form factor from BABAR
      !! Reference: Phys.Rev.D 86 (2012) 032013
      !! Input: squared momentum (q2)
      !! Output: complex form factor (ff), squared form factor (ff2)
      !! Energy range: 4*m_pi^2 < s < 9 GeV^2
      subroutine ffpi_babar(q2,ff,ff2)

      implicit none

      double precision q2,ff2
      double complex ff

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      double precision ampi
      common/hparam/ampi

      double precision m_n(0:3),g_n(0:3)
      double precision mod_n(0:3),arg_n(0:3)
      double complex f_n(0:3)

      double complex bw_const,bw_gs

      double precision s_max
      common/smax/s_max

      !! Maximum safe q^2 (GeV^2)
      s_max = 9.d0
      call ffpi_warning(q2,s_max)

      !! Charged pion mass (GeV)
      m_pi = ampi

      !! Omega mass and width (GeV)
      m_om = 0.78191d0
      g_om = 0.00813d0

      !! Masses of vector mesons (GeV)
      m_n(0) = 0.77502d0
      m_n(1) = 1.493d0
      m_n(2) = 1.861d0
      m_n(3) = 2.254d0

      !! Widths of vector mesons (GeV)
      g_n(0) = 0.14959d0
      g_n(1) = 0.427d0
      g_n(2) = 0.316d0
      g_n(3) = 0.109d0

      !! Resonance amplitudes

      mod_n(0) = 1.644d-3
      mod_n(1) = 0.158
      mod_n(2) = 0.068
      mod_n(3) = 0.0051

      arg_n(0) = -0.011d0
      arg_n(1) = 3.76
      arg_n(2) = 1.39
      arg_n(3) = 0.7

      f_n(0) = mod_n(0)*dcmplx(cos(arg_n(0)),sin(arg_n(0)))
      f_n(1) = mod_n(1)*dcmplx(cos(arg_n(1)),sin(arg_n(1)))
      f_n(2) = mod_n(2)*dcmplx(cos(arg_n(2)),sin(arg_n(2)))
      f_n(3) = mod_n(3)*dcmplx(cos(arg_n(3)),sin(arg_n(3)))

      !! Pion form factor
      ff = (bw_gs(m_n(0),g_n(0),q2) / (1.d0 + f_n(0))
     &   * (1.d0 + f_n(0) * bw_const(m_om,g_om,q2))
     &   + f_n(1) * bw_gs(m_n(1),g_n(1),q2)
     &   + f_n(2) * bw_gs(m_n(2),g_n(2),q2)
     &   + f_n(3) * bw_gs(m_n(3),g_n(3),q2))
     &   / (1 + f_n(1) + f_n(2) + f_n(3))

      !! Squared pion form factor
      ff2 = real(ff)**2 + aimag(ff)**2

      return
      end

      !!==============================================================!!

      !! Charged pion form factor from BES III
      !! Reference: Phys.Lett.B 753 (2016) 629–638
      !! Input: squared momentum (q2)
      !! Output: complex form factor (ff), squared form factor (ff2)
      !! Energy range: 4*m_pi^2 < s < 9 GeV^2
      subroutine ffpi_besiii(q2,ff,ff2)

      implicit none

      double precision q2,ff2
      double complex ff

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      double precision ampi
      common/hparam/ampi

      double precision m_n(0:3),g_n(0:3)
      double precision mod_n(0:3),arg_n(0:3)
      double complex f_n(0:3)

      double complex bw_const,bw_gs

      double precision s_max
      common/smax/s_max

      !! Maximum safe q^2 (GeV^2)
      s_max = 9.d0
      call ffpi_warning(q2,s_max)

      !! Charged pion mass (GeV)
      m_pi = ampi

      !! Omega mass and width (GeV)
      m_om = 0.78191d0
      g_om = 0.00813d0

      !! Masses of vector mesons (GeV)
      m_n(0) = 0.77658d0
      m_n(1) = 1.493d0
      m_n(2) = 1.861d0
      m_n(3) = 2.254d0

      !! Widths of vector mesons(GeV)
      g_n(0) = 0.15205d0
      g_n(1) = 0.427d0
      g_n(2) = 0.316d0
      g_n(3) = 0.109d0

      !! Resonance amplitudes

      mod_n(0) = 1.92d-3
      mod_n(1) = 0.158
      mod_n(2) = 0.068
      mod_n(3) = 0.0051

      arg_n(0) = 0.15d0
      arg_n(1) = 3.76
      arg_n(2) = 1.39
      arg_n(3) = 0.7

      f_n(0) = mod_n(0)*dcmplx(cos(arg_n(0)),sin(arg_n(0)))
      f_n(1) = mod_n(1)*dcmplx(cos(arg_n(1)),sin(arg_n(1)))
      f_n(2) = mod_n(2)*dcmplx(cos(arg_n(2)),sin(arg_n(2)))
      f_n(3) = mod_n(3)*dcmplx(cos(arg_n(3)),sin(arg_n(3)))

      !! Pion form factor
      ff = (bw_gs(m_n(0),g_n(0),q2) / (1.d0 + f_n(0))
     &   * (1.d0 + f_n(0) * bw_const(m_om,g_om,q2))
     &   + f_n(1) * bw_gs(m_n(1),g_n(1),q2)
     &   + f_n(2) * bw_gs(m_n(2),g_n(2),q2)
     &   + f_n(3) * bw_gs(m_n(3),g_n(3),q2))
     &   / (1 + f_n(1) + f_n(2) + f_n(3))

      !! Squared pion form factor
      ff2 = real(ff)**2 + aimag(ff)**2

      return
      end

      !!==============================================================!!

      !! Charged pion form factor from KLOE-2
      !! Reference: JHEP 03 (2018) 173
      !! Input: squared momentum (q2)
      !! Output: complex form factor (ff), squared form factor (ff2)
      !! Energy range: 4*m_pi^2 < s < 1 GeV^2
      !! Note: the KLOE-2 form factor has never been fitted by the
      !!       collaboration, this is just an indicative expression
      subroutine ffpi_kloe2(q2,ff,ff2)

      implicit none

      double precision q2,ff2
      double complex ff

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      double precision ampi
      common/hparam/ampi

      double precision m_n(0:1),g_n(0:1)
      double precision mod_n(0:1),arg_n(0:1)
      double complex f_n(0:1)

      double precision mod_om,arg_om,beta
      double complex del_om

      double complex bw_const,bw_gs

      double precision m_phi,g_phi
      double precision mod_phi,arg_phi
      double complex del_phi

      double precision s_max
      common/smax/s_max

      !! Maximum safe q^2 (GeV^2)
      s_max = 1.d0
      call ffpi_warning(q2,s_max)

      !! Charged pion mass (GeV)
      m_pi = ampi

      !! Omega mass and width (GeV)
      m_om = 0.7836351147d0
      g_om = 0.0095130784222d0

      !! Masses of vector mesons (GeV)
      m_n(0) = 0.7745790930242d0
      m_n(1) = 1.465d0

      !! Widths of vector mesons (GeV)
      g_n(0) = 0.1462368086199d0
      g_n(1) = 0.400d0

      !! Resonance amplitudes

      mod_n(1) = 0.08280137948d0
      arg_n(1) = 3.1415926535d0

      mod_om = 0.0015644121297d0
      arg_om = 0.1818523456d0

      f_n(1) = mod_n(1)*dcmplx(cos(arg_n(1)),sin(arg_n(1)))
      del_om = mod_om*dcmplx(cos(arg_om),sin(arg_om))

      !! Pion form factor
      ff = (bw_gs(m_n(0),g_n(0),q2)
     &   * (1 + del_om*q2/m_om**2 * bw_const(m_om,g_om,q2))
     &   + f_n(1) * bw_gs(m_n(1),g_n(1),q2)) / (1 + f_n(1))

      !! Squared pion form factor
      ff2 = real(ff)**2 + aimag(ff)**2

      return
      end

      !!==============================================================!!

      !! Charged pion form factor used for the Strong2020 report
      !! Based on experimental data before CMD-3
      !! Input: squared momentum (q2)
      !! Output: complex form factor (ff), squared form factor (ff2)
      !! Energy range: 4*m_pi^2 < s < 9 GeV^2
      subroutine ffpi_strong2020(q2,ff,ff2)

      implicit none

      double precision q2,ff2
      double complex ff

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      double precision ampi
      common/hparam/ampi

      double precision m_n(0:3),g_n(0:3)
      double precision mod_n(0:3),arg_n(0:3)
      double complex f_n(0:3)

      double precision mod_om,arg_om,beta
      double complex del_om

      double complex bw_const,bw_gs

      double precision m_phi,g_phi
      double precision mod_phi,arg_phi
      double complex del_phi

      double precision s_max
      common/smax/s_max

      !! Maximum safe q^2 (GeV^2)
      s_max = 9.d0
      call ffpi_warning(q2,s_max)

      !! Charged pion mass (GeV)
      m_pi = ampi

      !! Omega mass and width (GeV)
      m_om = 0.78248d0
      g_om = 0.00855d0

      !! Phi mass and width (GeV)
      m_phi = 1.01947d0
      g_phi = 0.00425d0

      !! Masses of vector mesons (GeV)
      m_n(0) = 0.77456d0
      m_n(1) = 1.4859d0
      m_n(2) = 1.8668d0
      m_n(3) = 2.2645d0

      !! Widths of vector mesons (GeV)
      g_n(0) = 0.14832d0
      g_n(1) = 0.37360d0
      g_n(2) = 0.30334d0
      g_n(3) = 0.11327d0

      !! Resonance amplitudes

      mod_n(1) = 0.14104d0
      mod_n(2) = 0.0614d0
      mod_n(3) = 0.0047d0
      mod_om = 0.00158d0
      mod_phi = 0.00045d0

      arg_n(1) = 3.7797d0
      arg_n(2) = 1.429d0
      arg_n(3) = 0.921d0
      arg_om = 0.075d0
      arg_phi = 2.888d0

      f_n(1) = mod_n(1)*dcmplx(cos(arg_n(1)),sin(arg_n(1)))
      f_n(2) = mod_n(2)*dcmplx(cos(arg_n(2)),sin(arg_n(2)))
      f_n(3) = mod_n(3)*dcmplx(cos(arg_n(3)),sin(arg_n(3)))
      del_phi = mod_phi*dcmplx(cos(arg_phi),sin(arg_phi))
      del_om = mod_om*dcmplx(cos(arg_om),sin(arg_om))

      !! Pion form factor
      ff = (bw_gs(m_n(0),g_n(0),q2)
     &   * (1 + del_om*q2/m_om**2 * bw_const(m_om,g_om,q2)
     &   + del_phi*q2/m_phi**2 * bw_const(m_phi,g_phi,q2))
     &   + f_n(1) * bw_gs(m_n(1),g_n(1),q2)
     &   + f_n(2) * bw_gs(m_n(2),g_n(2),q2)
     &   + f_n(3) * bw_gs(m_n(3),g_n(3),q2))
     &   / (1 + f_n(1) + f_n(2) + f_n(3))

!      ff=bw_const(0.77526d0,0.1474d0,q2)

      !! Squared pion form factor
      ff2 = real(ff)**2 + aimag(ff)**2

      return
      end

      !!==============================================================!!

      !! Charged pion form factor from dispersion relations
      !! Based on experimental data before CMD-3
      !! Reference: JHEP 02 (2019) 006
      !! Input: squared momentum (q2)
      !! Output: squared form factor (ff2)
      !! Energy range: 4*m_pi^2 < s < 4 GeV^2
      subroutine ffpi_bern(q2,ff,ff2)

      implicit none

      double precision q2,ff2
      double complex ff

      double complex ffpi_table,ffpi_interpolate
      character*50 input_file

      integer nitems,interpol

      !! Input file location
      input_file = 'ffpi_bern_input.txt'

      !! Number of items in the input table (< 2000)
      nitems = 1000

      !! Flag to activate the linear interpolation of the input table
      interpol = 1

      !! Pion form factor
      if (interpol.eq.1) then
         ff = ffpi_interpolate(sqrt(q2),input_file,nitems)
      else
         ff = ffpi_table(sqrt(q2),input_file,nitems)
      endif

      !! Squared pion form factor
      ff2 = real(ff)**2 + aimag(ff)**2

      return
      end

      !!==============================================================!!

      !! Pion form factor contribution from Dual QCD for n >= 5
      !! Reference: Phys.Lett.B 512 (2001) 331-334
      !! Input: momentum^2 (m_n), number of resonances (nmax)
      double complex function ff_dqcd(q2,nmax)

      implicit none

      double precision q2
      integer nmax

      double precision c_k,m_k,g_k
      double complex bw_gs
      integer k

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      integer nml
      parameter (nml = 3000)
      double precision c(nml),m(nml),g(nml)

      double precision beta,iloc
      common/par_dqcd/beta

      double precision pi
      integer ifirst
      common/ffqcdcmn/c,m,g,pi,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         if (nml.lt.nmax) then
            print*,'Error: Underdimensioned arrays in ff_dqcd function'
         endif
         pi = 4.d0*atan(1.d0)

         c = 0.d0
         m = 0.d0
         g = 0.d0

         !! Filling c(k), m(k), g(k)
         iloc = (6-beta)*(5-beta)*(4-beta)*(3-beta)*(2-beta)/120.d0
         do k=6,nmax
            iloc = iloc * (1.d0 - (beta-1.d0)/k)
            c_k = (-1.d0)**k * gamma(beta-0.5d0) * 2.d0/sqrt(pi)
     &           / (1.d0 + 2.d0*k) * iloc/pi * gamma(2.d0-beta)
     &           * sin(pi*(beta-1.d0-k))
            c(k) = c_k

            m_k = sqrt(m_rho**2 * (1.d0 + 2.d0*k))
            m(k) = m_k

            g_k = g_rho/m_rho*m_k
            g(k) = g_k
         enddo
         ifirst = 1
      endif

      ff_dqcd = 0.d0
      do k=6,nmax
         ff_dqcd = ff_dqcd + c(k) * bw_gs(m(k),g(k),q2)
      enddo

      return
      end

      !!==============================================================!!

      !! Breit-Wigner distribution by Gounaris-Sakurai
      !! Reference: Phys.Rev.Lett. 21 (1968) 244-247
      !! Input: hadron mass (m_n), decay width (g_n), momentum^2 (q2)
      double complex function bw_gs(m_n,g_n,q2)

      implicit none

      double precision m_n,g_n,q2

      double precision q,p_n,p_q
      double precision d_n,h_n,h_q,dh_n,f_n,g_n_q

      double precision m_pi,m_rho,g_rho,m_om,g_om
      common/par_hadrons/m_pi,m_rho,g_rho,m_om,g_om

      double precision pi
      double complex i

      pi = 4.d0*atan(1.d0)
      i = (0.d0,1.d0)

      q = sqrt(q2)
      p_n = sqrt(m_n**2 - 4.d0*m_pi**2) / 2.d0
      p_q = sqrt(q2 - 4.d0*m_pi**2) / 2.d0
! CARLO      p_q = sqrt(abs(q2 - 4.d0*m_pi**2)) / 2.d0

      d_n = 3.d0/pi*(m_pi/p_n)**2 * log((m_n + 2.d0*p_n)/(2.d0*m_pi))
     &    + m_n/(2.d0*pi*p_n) - m_n*m_pi**2/pi/p_n**3

      h_n = 2.d0/pi*p_n/m_n * log((m_n + 2d0*p_n)/(2.d0*m_pi))
      h_q = 2.d0/pi*p_q/q * log((q + 2d0*p_q)/(2.d0*m_pi))

      dh_n = h_n/8.d0*(1.d0/p_n**2 - 4/m_n**2) + 1/(2.d0*pi*m_n**2)

      g_n_q = g_n*m_n/q*(p_q/p_n)**3

      f_n = (h_q - h_n)*p_q**2 + (m_n**2 - q2)*dh_n*p_n**2
      f_n = f_n*g_n*m_n**2 / p_n**3

      bw_gs = (m_n**2 + d_n*m_n*g_n) / (m_n**2 - q2 + f_n - i*m_n*g_n_q)

      return
      end

      !!==============================================================!!

      !! Constant Breit-Wigner distribution
      !! Input: hadron mass (m_n), decay width (g_n), momentum^2 (q2)
      double complex function bw_const(m_n,g_n,q2)

      implicit none

      double precision m_n,g_n,q2

      double complex i
      i = (0.d0,1.d0)

      bw_const = m_n*m_n / (m_n*m_n - q2 - i*m_n*g_n)

      return
      end

      !!==============================================================!!

      !! Breit-Wigner distribution for GVMD
      !! Input: hadron mass (m_n), decay width (g_n), momentum^2 (q2)
      double complex function bw_gvmd(m_n,g_n,q2)

      implicit none

      double precision m_n,g_n,q2

      double complex i
      i = (0.d0,1.d0)

      bw_gvmd = (m_n*m_n - i*m_n*g_n) / (m_n*m_n - q2 - i*m_n*g_n)

      return
      end

      !!==============================================================!!

      !! Extract the squared form factor from an input table
      !! The value is selected according to the closest momentum
      !! Input: squared momentum (q2), input file name (file_name),
      !!        number of elements in the table (nitems)
      double complex function ffpi_table(sqs0,file_name,nitems)

      implicit none

      double precision sqs0
      character*50 file_name
      integer nitems

      !! Maximum number of items allowed in the input table
      integer nitems_max
      parameter (nitems_max = 2000)

      double precision sqs(1:nitems_max)
      double precision ffpi_re(1:nitems_max),ffpi_im(1:nitems_max)
      double precision ffpi_re_tmp(0:0),ffpi_im_tmp(0:0)
      integer i

      integer ifirst
      common/ffpitablecmn/sqs,ffpi_re,ffpi_im,ifirst
      data ifirst /0/

      if(nitems.gt.nitems_max) then
         print*, 'Error: Input table with too many items (nitems_max =',
     &           nitems_max, ')'
      endif

      !! The input table must have 3 columns:
      !!    1) Momentum (GeV)
      !!    2) Real part of the form factor
      !!    2) Imaginary part of the form factor
      if (ifirst.eq.0) then
         open(12,file=file_name,status='old')
         do i=1,nitems
            read(12,*) sqs(i),ffpi_re(i),ffpi_im(i)
         enddo
         close(12)
         ifirst = 1
      endif

      if (sqs0.lt.minval(sqs).or.sqs0.gt.maxval(sqs)) then
         print*, 'Error: Pion form factor out of energy range'
         stop
      endif

      !! Closest value
      ffpi_re_tmp = ffpi_re(minloc(abs(sqs-sqs0)))
      ffpi_im_tmp = ffpi_im(minloc(abs(sqs-sqs0)))

      ffpi_table = cmplx(ffpi_re_tmp(0),ffpi_im_tmp(0))

      return
      end

      !!==============================================================!!

      !! Extract the squared form factor from an input table
      !! The value is selected through linear interpolation
      !! Input: squared momentum (q2), input file name (file_name),
      !!        number of elements in the table (nitems)
      double complex function ffpi_interpolate(sqs0,file_name,nitems)

      implicit none

      double precision sqs0
      character*50 file_name
      integer nitems

      !! Maximum number of items allowed in the input table
      integer nitems_max
      parameter (nitems_max = 2000)

      double precision sqs(1:nitems_max),sqs_loc,sqs_diff
      double precision ffpi_re(1:nitems_max),ffpi_im(1:nitems_max)
      double complex ffpi_t(1:nitems_max)
      double complex ffpi_tmp(0:0),ffpi_tmp2(0:0)
      double precision sqs_tmp(0:0),sqs_tmp2(0:0)
      double precision x1,x2
      double complex y1,y2
      integer i

      integer ifirst
      common/ffpiintcmn/sqs,ffpi_re,ffpi_im,ffpi_t,ifirst
      data ifirst /0/

      double complex im
      im = (0.d0,1.d0)

      if(nitems.gt.nitems_max) then
         print*, 'Error: Input table with too many items (nitems_max =',
     &           nitems_max, ')'
      endif

      !! The input table must have 3 columns:
      !!    1) Momentum (GeV)
      !!    2) Real part of the form factor
      !!    2) Imaginary part of the form factor
      if (ifirst.eq.0) then
         open(12,file=file_name,status='old')
         do i=1,nitems
            read(12,*) sqs(i),ffpi_re(i),ffpi_im(i)
            ffpi_t(i) = ffpi_re(i) + im*ffpi_im(i)
         enddo
         close(12)
         ifirst = 1
      endif

      if (sqs0.lt.minval(sqs).or.sqs0.gt.maxval(sqs)) then
         print*, 'Error: Pion form factor out of energy range'
         stop
      endif

      !! Closest value
      ffpi_tmp = ffpi_t(minloc(abs(sqs-sqs0)))
      sqs_tmp = sqs(minloc(abs(sqs-sqs0)))

      sqs_loc = sqs_tmp(0)
      sqs_diff = sqs0 - sqs_loc

      x1 = sqs_loc
      y1 = ffpi_tmp(0)

      !! Linear interpolation
      if(sqs_diff.gt.0) then
         sqs_tmp2 = sqs(minloc(abs(sqs-sqs0))+1)
         x2 = sqs_tmp2(0)
         ffpi_tmp2 = ffpi_t(minloc(abs(sqs-sqs0))+1)
         y2 = ffpi_tmp2(0)
         ffpi_interpolate = y1 + (sqs0-x1)*(y2-y1)/(x2-x1)
      else
         sqs_tmp2 = sqs(minloc(abs(sqs-sqs0))-1)
         x2 = sqs_tmp2(0)
         ffpi_tmp2 = ffpi_t(minloc(abs(sqs-sqs0))-1)
         y2 = ffpi_tmp2(0)
         ffpi_interpolate = y1 + (sqs0-x1)*(y2-y1)/(x2-x1)
      endif

      return
      end

      !!==============================================================!!

      !! Verify the form factor sum rule for the dispersive approach
      !! Input: maximum squared momentum in GeV^2 (q2max)
      !! Output: dispersive sum rule integral (sum)
      subroutine sum_rule(q2max,sum)

      implicit none

      double precision q2min,q2max,sum

      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi

      double precision ampi
      common/hparam/ampi

      double complex ffx
      double precision x,ffx2,func
      integer j,nstep

      double precision pi
      double complex i

      pi = 4.d0*atan(1.d0)
      i = (0.d0,1.d0)

      !! Numerical integral with the Simpson rule from 4*mpi^2 to q2max
      nstep = int(100*q2max)
      sum = 0.d0
      q2min = 4.d0*ampi**2 + 1d-12 ! Add 1 keV^2 for safety

      do j=0,nstep

         x = q2min + j*(q2max - q2min)/nstep
         call ffpi(x,ffx,ffx2)
         func = aimag(ffx)/x/pi

         if (j.eq.0) then
            sum = sum + (q2max-q2min)/(3.d0*nstep)*func
         elseif (mod(j,2).ne.0) then
            sum = sum + 4.d0*(q2max-q2min)/(3.d0*nstep)*func
         elseif (mod(j,2).eq.0.and.j.ne.0) then
            sum = sum + 2.d0*(q2max-q2min)/(3.d0*nstep)*func
         endif

      enddo

      return
      end

      !!==============================================================!!

      !! Subroutine to manage warnings on the pion form factor range
      !! Input: squared momentum (q2), max squared momentum (q2_max)
      subroutine ffpi_warning(q2,q2_max)

      implicit none

      double precision q2,q2_max

      integer iffpiwarning,warn_max
      data iffpiwarning /0/

      warn_max = 10

      if(q2.gt.q2_max.and.iffpiwarning.lt.warn_max) then
         print*, 'Warning: Pion form factor out of energy range (s =',
     &           q2, 'GeV^2 > s_max =', q2_max, 'GeV^2)'

         iffpiwarning = iffpiwarning + 1

         if(iffpiwarning.eq.warn_max) then
            print*, 'Suppressing further warnings on the pion form ',
     &              'factor range... be careful!'
         endif

      endif

      return
      end

      !!==============================================================!!
