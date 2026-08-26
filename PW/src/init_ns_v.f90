SUBROUTINE init_ns_v

   USE kinds,      ONLY : DP
   USE ions_base,  ONLY : nat, ityp, ntyp => nsp
   USE lsda_mod,   ONLY : nspin, starting_magnetization
   USE ldaV,       ONLY : ehub_lm_max, ehub_lm, ehub_l, &
                          ehub_l_choice, lacbn0, acbn0_type
   USE scf,        ONLY : rho
   USE uspp_param, ONLY : upf

   IMPLICIT NONE

   REAL(DP) :: totocc(3,ntyp)

   LOGICAL  :: nm

   INTEGER  :: na, nt
   INTEGER  :: nw, l, m, ldim
   INTEGER  :: is, majs, mins

   totocc(:,:) = 0.0_DP
   !
   IF ( lacbn0 .AND. acbn0_type == 1 ) THEN
      rho%na0(:,:,:,:) = 0.0_DP
   ELSE
      rho%ne0(:,:,:,:,:,:) = 0.0_DP
   ENDIF
   !
   DO nt = 1, ntyp
      !
      IF ( ehub_lm(nt) == 0 ) CYCLE
      !
      DO nw = 1, upf(nt)%nwfc
         !
         IF ( upf(nt)%oc(nw) >= 0.0_DP .AND. ehub_l_choice(nt,nw) == 1 ) THEN
            l = upf(nt)%lchi(nw)
            totocc(l+1,nt) = upf(nt)%oc(nw)
         ENDIF
         !
      ENDDO ! nw
      !
   ENDDO ! nt
   !
   DO na = 1, nat
      nt = ityp(na)
      !
      IF ( ehub_lm(nt) == 0 ) CYCLE
      !
      nm = .true.
      !
      IF ( nspin == 2 ) THEN
         IF ( starting_magnetization(nt) > 0.0_DP ) THEN
            nm = .false.
            majs = 1
            mins = 2
         ELSE IF ( starting_magnetization(nt) < 0.0_DP ) THEN
            nm = .false.
            majs = 2
            mins = 1
         ENDIF
      ENDIF
      !
      IF ( .NOT. nm ) THEN
         !
         DO m = 1, ehub_lm(nt)
            ldim = 2 * (ehub_l(m,nt) - 1) + 1
            !
            IF ( totocc(ehub_l(m,nt),nt) > ldim ) THEN
               !
               IF ( lacbn0 .AND. acbn0_type == 1 ) THEN
                  rho%na0(m,m,majs,na) = 1.0_DP
                  rho%na0(m,m,mins,na) = (totocc(ehub_l(m,nt),nt) - ldim) / ldim
               ELSE
                  rho%ne0(m,m,majs,na,na,1) = 1.0_DP
                  rho%ne0(m,m,mins,na,na,1) = (totocc(ehub_l(m,nt),nt) - ldim) / ldim
               ENDIF
            ELSE
               IF ( lacbn0 .AND. acbn0_type == 1 ) THEN
                  rho%na0(m,m,majs,na) = totocc(ehub_l(m,nt),nt) / ldim
               ELSE
                  rho%ne0(m,m,majs,na,na,1) = totocc(ehub_l(m,nt),nt) / ldim
               ENDIF
            ENDIF
            !
         ENDDO
         !
      ELSE
         !
         DO is = 1, nspin
            DO m = 1, ehub_lm(nt)
               ldim = 2 * (ehub_l(m,nt) - 1) + 1
               !
               IF ( lacbn0 .AND. acbn0_type == 1 ) THEN
                  rho%na0(m,m,is,na) = totocc(ehub_l(m,nt),nt) / 2.0_DP / ldim
               ELSE
                  rho%ne0(m,m,is,na,na,1) = totocc(ehub_l(m,nt),nt) / 2.0_DP / ldim
               ENDIF
               !
            ENDDO
         ENDDO
         !
      ENDIF
      !
   ENDDO
   !
   RETURN

END SUBROUTINE init_ns_v
!
! WIY
SUBROUTINE init_ns_v_nc
   !
   !! Noncollinear version of init_ns (WIY)
   !
   USE kinds,              ONLY : DP
   USE ions_base,          ONLY : nat, ityp, ntyp => nsp
   USE lsda_mod,           ONLY : nspin, starting_magnetization
   USE ldaV,               ONLY : ehub_lm_max, ehub_lm, ehub_l, &
                                  ehub_l_choice, lacbn0, acbn0_type
   USE noncollin_module,   ONLY : angle1, angle2
   USE scf,                ONLY : rho
   USE uspp_param,         ONLY : upf
   !
   IMPLICIT NONE
   !
   REAL(DP) :: totocc(3,ntyp), cosin
   COMPLEX(DP) :: esin, ns(4)
   !! parameter for rotating occp. matrix
   COMPLEX(DP) :: n_charge
   !! electron charge
   COMPLEX(DP) :: m_moment
   !! magnetic moment
   !
   LOGICAL  :: nm 
   !! true if the atom is nonmagnetic
   !
   INTEGER  :: na, nt
   INTEGER  :: nw, l, m, ldim
   INTEGER  :: is, majs, mins
   !
   totocc(:,:) = 0.0_DP
   !
   IF ( lacbn0 .AND. acbn0_type == 1 ) THEN
      rho%na0_nc(:,:,:,:) = 0.0_DP
   ELSE
      rho%ne0_nc(:,:,:,:,:,:) = 0.0_DP
   ENDIF
   !
   DO nt = 1, ntyp
      !
      IF ( ehub_lm(nt) == 0 ) CYCLE
      !
      DO nw = 1, upf(nt)%nwfc
         !
         l = upf(nt)%lchi(nw)
         !
         IF ( upf(nt)%oc(nw) >= 0.0_DP .AND. ehub_l_choice(nt,l+1) == nw ) THEN
            !
            totocc(l+1,nt) = upf(nt)%oc(nw)
            !
            IF ( (nw-1 > 0) .AND. ( trim(upf(nt)%els(nw)) == trim(upf(nt)%els(nw-1)) ))  THEN
            !
            totocc(l+1,nt) = totocc(l+1,nt) + upf(nt)%oc(nw-1)
            !
            ENDIF
            !
         ENDIF
         !
      ENDDO ! nw
      !
   ENDDO ! nt
   !
   DO na = 1, nat
      nt = ityp(na)
      !
      IF ( ehub_lm(nt) == 0 ) CYCLE
      !
      nm = .true.
      !
      IF ( starting_magnetization(nt) > 0.0_DP ) THEN
         nm = .false.
         majs = 1
         mins = 2
      ELSE IF ( starting_magnetization(nt) < 0.0_DP ) THEN
         nm = .false.
         majs = 2
         mins = 1
      ENDIF
      !
      IF ( .NOT. nm ) THEN
         !
         cosin = COS(angle1(nt))
         esin  = ( COS(angle2(nt)) + (0.0_DP,1.0_DP)*SIN(angle2(nt))) * SIN(angle1(nt))
         !
         DO m = 1, ehub_lm(nt)
            ldim = 2 * (ehub_l(m,nt) - 1) + 1
            !
            ! occp. matrix in quantized axis
            IF ( totocc(ehub_l(m,nt),nt) > ldim ) THEN
               ns(majs) = 1.0_DP 
               ns(mins) = (totocc(ehub_l(m,nt),nt) - ldim) / ldim
            ELSE
               ns(majs) = totocc(ehub_l(m,nt),nt) / ldim 
               ns(mins) = 0.0_DP 
            ENDIF
            !
            n_charge = ns(1) + ns(2)
            m_moment = ns(1) - ns(2) 
            !
            ! rotating occp. matrix
            ns(1) = ( n_charge + m_moment * cosin) / 2.0_DP
            ns(2) = m_moment * esin / 2.0_DP 
            ns(3) = m_moment * CONJG (esin) / 2.0_DP 
            ns(4) = ( n_charge - m_moment * cosin) / 2.0_DP 
            !
            IF ( lacbn0 .AND. acbn0_type == 1 ) THEN
               rho%na0_nc(m,m,:,na) = ns(:)
            ELSE
               rho%ne0_nc(m,m,:,na,na,1) = ns(:)
            ENDIF
            !
         ENDDO
         !
      ELSE
         !
         DO m = 1, ehub_lm(nt)
            ldim = 2 * (ehub_l(m,nt) - 1) + 1
            !
            IF ( lacbn0 .AND. acbn0_type == 1 ) THEN
               rho%na0_nc(m,m,1,na) = totocc(ehub_l(m,nt),nt) / 2.0_DP / ldim
               rho%na0_nc(m,m,4,na) = totocc(ehub_l(m,nt),nt) / 2.0_DP / ldim              
            ELSE
               rho%ne0_nc(m,m,1,na,na,1) = totocc(ehub_l(m,nt),nt) / 2.0_DP / ldim
               rho%ne0_nc(m,m,4,na,na,1) = totocc(ehub_l(m,nt),nt) / 2.0_DP / ldim
            ENDIF
            !
         ENDDO
         !
      ENDIF
      !
   ENDDO
   !
   RETURN

END SUBROUTINE init_ns_v_nc
! WIY
